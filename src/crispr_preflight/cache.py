"""Polite rate-limited HTTP client with local SQLite TTL caching and gzip compression.

Respects fair-use policies across public biological APIs (Ensembl <=15 req/s,
gnomAD <=10 req/min, Open Targets, GTEx, and Human Protein Atlas) with
automatic HTTP 429 Retry-After handling, gzip decompression, and exponential backoff.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import sqlite3
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Optional

USER_AGENT = "crispr-preflight/0.1.0 (Open-Source Wet-Lab CRISPR Flight Check; Python urllib)"
DEFAULT_TTL_SECONDS = 7 * 24 * 3600  # 7 days


class NetworkFetchError(RuntimeError):
    """Raised when an HTTP endpoint fails after retries."""


def default_cache_path() -> Path:
    """Return platform-appropriate cache SQLite file path."""
    override = os.environ.get("CRISPR_PREFLIGHT_CACHE")
    if override:
        return Path(override)
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else (Path.home() / ".cache")
    return base / "crispr-preflight" / "http_cache_v1.sqlite3"


class CachedHttpClient:
    """Rate-limited HTTP/GraphQL client backed by a local SQLite TTL cache."""

    def __init__(
        self,
        cache_path: Optional[Path] = None,
        ttl_seconds: int = DEFAULT_TTL_SECONDS,
        use_cache: bool = True,
        timeout_seconds: float = 14.0,
    ) -> None:
        self.use_cache = use_cache
        self.ttl_seconds = ttl_seconds
        self.timeout_seconds = timeout_seconds
        self.cache_path = cache_path or default_cache_path()
        self._db_lock = threading.Lock()
        self._host_locks_guard = threading.Lock()
        self._host_locks: dict[str, threading.Lock] = {}
        self._last_request_by_host: dict[str, float] = {}
        self._mem_cache: dict[str, tuple[bytes, float]] = {}
        self._min_interval_by_host: dict[str, float] = {
            "rest.ensembl.org": 0.07,              # <= 14.2 req/s (limit is 15/s)
            "www.ebi.ac.uk": 0.08,                 # EMBL-EBI Proteins & InterPro APIs
            "rest.uniprot.org": 0.08,              # UniProt REST API
            "eutils.ncbi.nlm.nih.gov": 0.34,       # NCBI E-utilities (<= 3 req/s without key)
            "gnomad.broadinstitute.org": 6.1,      # <= 10 req/min
            "api.platform.opentargets.org": 0.12,  # GraphQL v4 handles burst queries cleanly
            "gtexportal.org": 0.15,                # <= ~6.6 req/s
            "www.proteinatlas.org": 0.20,          # <= 5 req/s
            "www.cbioportal.org": 0.08,            # cBioPortal CCLE REST API
            "mygene.info": 0.08,                   # BioThings MyGene.info v3
            "myvariant.info": 0.08,                # BioThings MyVariant.info v1
        }
        if self.use_cache:
            self._init_db()

    def _init_db(self) -> None:
        create_sql = """
            CREATE TABLE IF NOT EXISTS responses (
                cache_key TEXT PRIMARY KEY,
                url TEXT NOT NULL,
                status_code INTEGER NOT NULL,
                body BLOB NOT NULL,
                fetched_at REAL NOT NULL
            )
        """
        for candidate in (
            self.cache_path,
            Path(tempfile.gettempdir()) / "crispr-preflight" / "http_cache_v1.sqlite3",
        ):
            try:
                candidate.parent.mkdir(parents=True, exist_ok=True)
                with sqlite3.connect(candidate) as conn:
                    conn.execute(create_sql)
                    conn.commit()
                self.cache_path = candidate
                return
            except (OSError, sqlite3.Error):
                continue
        # Read-only container fallback: keep in-memory cache active without failing startup
        self.use_cache = False

    @staticmethod
    def _make_key(method: str, url: str, headers: dict[str, str], body: Optional[bytes]) -> str:
        accept = headers.get("Accept", "")
        ctype = headers.get("Content-Type", "")
        raw = f"{method.upper()}|{url}|{accept}|{ctype}|".encode("utf-8") + (body or b"")
        return hashlib.sha256(raw).hexdigest()

    def _get_cached(self, key: str) -> Optional[bytes]:
        if not self.use_cache:
            return None
        now = time.time()
        mem_hit = self._mem_cache.get(key)
        if mem_hit is not None:
            body, fetched_at = mem_hit
            if now - fetched_at <= self.ttl_seconds:
                return body
        if not self.cache_path.exists():
            return None
        with self._db_lock, sqlite3.connect(self.cache_path) as conn:
            row = conn.execute(
                "SELECT body, fetched_at FROM responses WHERE cache_key = ?",
                (key,),
            ).fetchone()
        if row is None:
            return None
        body, fetched_at = row
        if now - float(fetched_at) > self.ttl_seconds:
            return None
        raw_bytes = bytes(body)
        if len(self._mem_cache) < 1024:
            self._mem_cache[key] = (raw_bytes, float(fetched_at))
        return raw_bytes

    def _put_cached(self, key: str, url: str, status_code: int, body: bytes) -> None:
        if not self.use_cache:
            return
        now = time.time()
        if len(self._mem_cache) >= 1024:
            self._mem_cache.clear()
        self._mem_cache[key] = (body, now)
        with self._db_lock, sqlite3.connect(self.cache_path) as conn:
            conn.execute(
                """
                INSERT OR REPLACE INTO responses (cache_key, url, status_code, body, fetched_at)
                VALUES (?, ?, ?, ?, ?)
                """,
                (key, url, status_code, body, now),
            )
            conn.commit()

    def clear_cache(self) -> int:
        """Delete all cached entries and return number of removed rows."""
        self._mem_cache.clear()
        if not self.cache_path.exists():
            return 0
        with self._db_lock, sqlite3.connect(self.cache_path) as conn:
            cur = conn.execute("DELETE FROM responses")
            conn.commit()
            return cur.rowcount or 0

    def cache_stats(self) -> dict[str, Any]:
        """Return cache file path and entry count."""
        if not self.cache_path.exists():
            return {"path": str(self.cache_path), "entries": 0, "exists": False}
        with self._db_lock, sqlite3.connect(self.cache_path) as conn:
            count = conn.execute("SELECT COUNT(*) FROM responses").fetchone()[0]
        return {"path": str(self.cache_path), "entries": int(count), "exists": True}

    def _throttle(self, url: str) -> None:
        host = urllib.parse.urlparse(url).netloc.lower()
        with self._host_locks_guard:
            hlock = self._host_locks.setdefault(host, threading.Lock())
        min_interval = self._min_interval_by_host.get(host, 0.10)
        with hlock:
            last = self._last_request_by_host.get(host)
            if last is not None:
                elapsed = time.monotonic() - last
                if elapsed < min_interval:
                    time.sleep(min_interval - elapsed)
            self._last_request_by_host[host] = time.monotonic()

    def fetch_bytes(
        self,
        url: str,
        *,
        method: str = "GET",
        headers: Optional[dict[str, str]] = None,
        json_body: Optional[dict[str, Any]] = None,
        max_retries: int = 2,
        timeout_override: Optional[float] = None,
    ) -> bytes:
        """Fetch raw response bytes with gzip compression, caching, rate limiting, and backoff."""
        req_headers = {
            "User-Agent": USER_AGENT,
            "Accept-Encoding": "gzip",
        }
        if headers:
            req_headers.update(headers)
        encoded_body: Optional[bytes] = None
        if json_body is not None:
            encoded_body = json.dumps(json_body, sort_keys=True).encode("utf-8")
            req_headers.setdefault("Content-Type", "application/json")

        key = self._make_key(method, url, req_headers, encoded_body)
        cached = self._get_cached(key)
        if cached is not None:
            return cached

        eff_timeout = timeout_override if timeout_override is not None else self.timeout_seconds
        backoff = 1.0
        last_exc: Optional[Exception] = None
        for attempt in range(max_retries + 1):
            self._throttle(url)
            req = urllib.request.Request(
                url,
                data=encoded_body,
                headers=req_headers,
                method=method.upper(),
            )
            try:
                with urllib.request.urlopen(req, timeout=eff_timeout) as resp:
                    status = getattr(resp, "status", 200)
                    payload = resp.read()
                    content_enc = (resp.headers.get("Content-Encoding") or "").lower()
                    if ("gzip" in content_enc or payload[:2] == b"\x1f\x8b") and payload:
                        payload = gzip.decompress(payload)
                    if 200 <= status < 300:
                        self._put_cached(key, url, status, payload)
                    return payload
            except urllib.error.HTTPError as exc:
                last_exc = exc
                if exc.code == 429 or 500 <= exc.code < 600:
                    if attempt < max_retries:
                        retry_after = exc.headers.get("Retry-After") if exc.headers else None
                        wait_s = float(retry_after) if (retry_after and retry_after.isdigit()) else backoff
                        time.sleep(min(wait_s, 12.0))
                        backoff *= 2.0
                        continue
                raise
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                last_exc = exc
                if attempt < max_retries:
                    time.sleep(backoff)
                    backoff *= 2.0
                    continue
                raise NetworkFetchError(f"Network timeout/error fetching {url}: {exc}") from exc

        if last_exc is not None:
            raise NetworkFetchError(f"Failed fetching {url}: {last_exc}") from last_exc
        raise NetworkFetchError(f"Failed to fetch {url}")

    def fetch_json(
        self,
        url: str,
        *,
        method: str = "GET",
        headers: Optional[dict[str, str]] = None,
        json_body: Optional[dict[str, Any]] = None,
        max_retries: int = 2,
        timeout_override: Optional[float] = None,
    ) -> Any:
        """Fetch and decode a JSON endpoint."""
        req_headers = {"Accept": "application/json"}
        if headers:
            req_headers.update(headers)
        raw = self.fetch_bytes(
            url,
            method=method,
            headers=req_headers,
            json_body=json_body,
            max_retries=max_retries,
            timeout_override=timeout_override,
        )
        if not raw:
            return None
        return json.loads(raw.decode("utf-8"))

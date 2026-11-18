"""Unit tests for the crispr-preflight web application HTML & API routes."""

from __future__ import annotations

import io
import json
import re
import shutil
import subprocess
import tempfile
import unittest

from unittest.mock import patch

from crispr_preflight.web import WEB_APP_HTML, make_handler


class _DummySocket:
    def __init__(self, request_bytes: bytes) -> None:
        self._rfile = io.BytesIO(request_bytes)
        self._wfile = io.BytesIO()

    def makefile(self, mode: str, *args: object) -> io.BytesIO:
        return self._rfile if "r" in mode else self._wfile

    def sendall(self, data: bytes) -> None:
        self._wfile.write(data)


class TestWebAppRoutes(unittest.TestCase):
    def _dispatch_get(self, path: str) -> tuple[str, bytes]:
        handler_cls = make_handler(use_cache=True)
        raw_req = f"GET {path} HTTP/1.1\r\nHost: localhost\r\n\r\n".encode("utf-8")
        sock = _DummySocket(raw_req)
        handler_cls(sock, ("127.0.0.1", 12345), None)  # type: ignore[arg-type]
        raw_resp = sock._wfile.getvalue()
        head, _, body = raw_resp.partition(b"\r\n\r\n")
        return head.decode("utf-8"), body

    def test_index_html_served(self) -> None:
        head, body = self._dispatch_get("/")
        self.assertIn("200 OK", head)
        self.assertIn("Run Flight Check", body.decode("utf-8"))
        self.assertEqual(len(body), len(WEB_APP_HTML.encode("utf-8")))

    def test_health_endpoint(self) -> None:
        head, body = self._dispatch_get("/api/health")
        self.assertIn("200 OK", head)
        data = json.loads(body.decode("utf-8"))
        self.assertEqual(data["status"], "ok")
        self.assertNotIn("path", data.get("cache", {}))
        self.assertNotIn("/Users/", body.decode("utf-8"))

    def test_internal_os_error_never_leaks_local_filepath(self) -> None:
        with patch(
            "crispr_preflight.engine.runner.PreflightRunner.run",
            side_effect=PermissionError(13, "Permission denied", "/Users/private_user/.cache/db.sqlite3"),
        ):
            head, body = self._dispatch_get("/api/check?gene=STAT3")
        text = body.decode("utf-8")
        self.assertIn("500", head)
        self.assertNotIn("private_user", text)
        self.assertNotIn("/Users/", text)
        self.assertNotIn(".sqlite3", text)

    def test_missing_gene_returns_400(self) -> None:
        head, body = self._dispatch_get("/api/check")
        self.assertIn("400", head)
        data = json.loads(body.decode("utf-8"))
        self.assertIn("error", data)

    def test_invalid_species_returns_400(self) -> None:
        head, body = self._dispatch_get("/api/check?species=martian&gene=Trp53")
        self.assertIn("400", head)
        data = json.loads(body.decode("utf-8"))
        self.assertIn("Unsupported species", data.get("error", ""))

    def test_invalid_nuclease_returns_400(self) -> None:
        head, body = self._dispatch_get("/api/check?species=mouse&nuclease=talen99&gene=Trp53")
        self.assertIn("400", head)
        data = json.loads(body.decode("utf-8"))
        self.assertIn("Unsupported nuclease", data.get("error", ""))

    def test_rate_limiter_and_uri_length_guard(self) -> None:
        # 1. Oversized URI (>2048 chars) returns 414 immediately
        long_path = "/api/check?gene=" + ("A" * 2200)
        head_long, _ = self._dispatch_get(long_path)
        self.assertIn("414", head_long)

        # 2. Sliding-window per-client rate limiter returns 429 + Retry-After when exceeded
        handler_cls = make_handler(use_cache=False, rate_limit_per_min=2)
        statuses: list[str] = []
        for i in range(3):
            raw_req = f"GET /api/snp-check?specs=%5B%5D&n={i} HTTP/1.1\r\nHost: localhost\r\n\r\n".encode("utf-8")
            sock = _DummySocket(raw_req)
            handler_cls(sock, ("127.0.0.1", 12345), None)  # type: ignore[arg-type]
            h, _, _ = sock._wfile.getvalue().partition(b"\r\n\r\n")
            statuses.append(h.decode("utf-8"))
        self.assertIn("200 OK", statuses[0])
        self.assertIn("200 OK", statuses[1])
        self.assertIn("429", statuses[2])
        self.assertIn("Retry-After: 30", statuses[2])

    def test_vercel_serverless_entrypoint_routing(self) -> None:
        from api.index import handler as vercel_handler

        raw_req = b"GET /api/index?__route=health HTTP/1.1\r\nHost: crispr-preflight.vercel.app\r\n\r\n"
        sock = _DummySocket(raw_req)
        vercel_handler(sock, ("127.0.0.1", 12345), None)  # type: ignore[arg-type]
        head, _, body = sock._wfile.getvalue().partition(b"\r\n\r\n")
        self.assertIn("200 OK", head.decode("utf-8"))
        self.assertEqual(json.loads(body.decode("utf-8"))["status"], "ok")

    def test_snp_check_empty_specs_returns_200(self) -> None:
        head, body = self._dispatch_get("/api/snp-check?specs=%5B%5D")
        self.assertIn("200 OK", head)
        data = json.loads(body.decode("utf-8"))
        self.assertEqual(data, {"results": {}})

        # Also test POST /api/snp-check with empty specs (zero network calls)
        handler_cls = make_handler(use_cache=True)
        payload = b'{"specs": []}'
        raw_post = (
            b"POST /api/snp-check HTTP/1.1\r\n"
            b"Host: localhost\r\n"
            b"Content-Type: application/json\r\n"
            + f"Content-Length: {len(payload)}\r\n\r\n".encode("utf-8")
            + payload
        )
        sock = _DummySocket(raw_post)
        handler_cls(sock, ("127.0.0.1", 12345), None)  # type: ignore[arg-type]
        raw_resp = sock._wfile.getvalue()
        p_head, _, p_body = raw_resp.partition(b"\r\n\r\n")
        self.assertIn("200 OK", p_head.decode("utf-8"))
        self.assertEqual(json.loads(p_body.decode("utf-8")), {"results": {}})

    def test_embedded_javascript_has_valid_syntax(self) -> None:
        """Ensure the <script> block inside WEB_APP_HTML parses cleanly with zero JS SyntaxErrors."""
        m = re.search(r"<script>([\s\S]*?)</script>", WEB_APP_HTML)
        self.assertIsNotNone(m)
        assert m is not None
        js_code = m.group(1)
        # Guard against unescaped Python literal newlines inside single-quoted JS strings
        self.assertNotIn(".join('\n')", js_code)

        node_bin = shutil.which("node")
        if node_bin:
            with tempfile.NamedTemporaryFile("w", suffix=".js", encoding="utf-8") as tf:
                tf.write(js_code)
                tf.flush()
                proc = subprocess.run([node_bin, "--check", tf.name], capture_output=True, text=True)
                self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        elif shutil.which("osascript"):
            # macOS JavaScriptCore syntax check via Function constructor
            check_expr = f"new Function({json.dumps(js_code)}); 'OK';"
            proc = subprocess.run(
                ["osascript", "-l", "JavaScript", "-e", check_expr],
                capture_output=True,
                text=True,
            )
            self.assertEqual(proc.returncode, 0, msg=proc.stderr)


if __name__ == "__main__":
    unittest.main()

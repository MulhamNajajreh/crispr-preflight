"""Vercel serverless Python entrypoint routing `/`, `/api/check`, `/api/snp-check`, and `/api/health`."""

from __future__ import annotations

import sys
import urllib.parse
from pathlib import Path

_SRC = str(Path(__file__).resolve().parent.parent / "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from crispr_preflight.web import make_handler  # noqa: E402

_BaseHandler = make_handler(use_cache=True)


class handler(_BaseHandler):  # noqa: N801
    def do_GET(self) -> None:  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        qs = urllib.parse.parse_qs(parsed.query)
        route_hint = (qs.get("_route") or qs.get("__route") or [""])[0].strip()
        if route_hint in ("check", "snp-check", "health"):
            self.path = f"/api/{route_hint}?{parsed.query}" if parsed.query else f"/api/{route_hint}"
        elif parsed.path in ("", "/api", "/api/index", "/api/index.py"):
            self.path = f"/?{parsed.query}" if parsed.query else "/"
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        parsed = urllib.parse.urlparse(self.path)
        if parsed.path in ("", "/", "/api", "/api/index", "/api/index.py"):
            self.path = f"/api/snp-check?{parsed.query}" if parsed.query else "/api/snp-check"
        super().do_POST()

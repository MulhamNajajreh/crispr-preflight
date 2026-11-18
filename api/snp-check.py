"""Vercel serverless Python entrypoint for `/api/snp-check`."""

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
        qs = urllib.parse.urlparse(self.path).query
        self.path = f"/api/snp-check?{qs}" if qs else "/api/snp-check"
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        self.path = "/api/snp-check"
        super().do_POST()

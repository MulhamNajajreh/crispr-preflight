"""Vercel serverless Python entrypoint for `/api/health`."""

from __future__ import annotations

import sys
from pathlib import Path

_SRC = str(Path(__file__).resolve().parent.parent / "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from crispr_preflight.web import make_handler  # noqa: E402

_BaseHandler = make_handler(use_cache=True)


class handler(_BaseHandler):  # noqa: N801
    def do_GET(self) -> None:  # noqa: N802
        self.path = "/api/health"
        super().do_GET()

"""Allow running via `python3 -m crispr_preflight`."""

from crispr_preflight.cli import main

if __name__ == "__main__":
    raise SystemExit(main())

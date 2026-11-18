"""Command-line interface for crispr-preflight (`check`, `batch`, `serve`, and `cache`)."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Optional

from crispr_preflight import __version__
from crispr_preflight.cache import CachedHttpClient
from crispr_preflight.engine.runner import PreflightRunner
from crispr_preflight.renderers.html_report import render_html
from crispr_preflight.renderers.markdown import render_markdown
from crispr_preflight.renderers.terminal import render_terminal
from crispr_preflight.web import run_web_server


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="crispr-preflight",
        description=(
            "Pre-experiment biological flight check for CRISPR knockout, CRISPRi, "
            "and target validation experiments."
        ),
    )
    parser.add_argument("--version", action="version", version=f"crispr-preflight {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- check ---
    p_check = subparsers.add_parser(
        "check",
        help="Run a biological preflight check for a single human gene",
    )
    p_check.add_argument("gene", help="HGNC gene symbol, alias (e.g. BRG1), or Ensembl ID (ENSG...)")
    p_check.add_argument(
        "-t", "--tissue", default=None, help="Target tissue/lineage context (e.g. lung, blood, brain, colon)"
    )
    p_check.add_argument(
        "-c", "--cell-line", default=None, help="DepMap cell line name (e.g. A549, K562, MCF7, HeLa)"
    )
    p_check.add_argument("--all-exons", action="store_true", help="Display all coding exons in 5'->3' order")
    p_check.add_argument("--no-pam", action="store_true", help="Skip SpCas9 NGG PAM sequence scan")
    p_check.add_argument("--no-cache", action="store_true", help="Bypass local SQLite HTTP cache")
    p_check.add_argument("--html", metavar="PATH", help="Save standalone interactive HTML+SVG report")
    p_check.add_argument("--markdown", metavar="PATH", help="Save Markdown lab-notebook dossier")
    p_check.add_argument("--json", metavar="PATH", help="Save full structured JSON output")
    p_check.add_argument("--quiet", action="store_true", help="Suppress terminal report when writing files")

    # --- batch ---
    p_batch = subparsers.add_parser(
        "batch",
        help="Audit a list of target genes (comma-separated or newline-delimited text file)",
    )
    p_batch.add_argument(
        "targets",
        help="Path to text file with one gene per line, or comma-separated symbols (e.g. SMARCA4,SOX2,PLK1)",
    )
    p_batch.add_argument("-t", "--tissue", default=None, help="Target tissue/lineage context")
    p_batch.add_argument("-c", "--cell-line", default=None, help="DepMap cell line name")
    p_batch.add_argument(
        "-o", "--output", default="preflight_batch_summary.csv", help="Output CSV summary table path"
    )
    p_batch.add_argument("--html-dir", default=None, help="Optional directory to write per-gene HTML reports")
    p_batch.add_argument("--no-pam", action="store_true", help="Skip SpCas9 NGG PAM sequence scan")
    p_batch.add_argument("--no-cache", action="store_true", help="Bypass local SQLite HTTP cache")

    # --- cache ---
    p_cache = subparsers.add_parser("cache", help="Inspect or clear the local SQLite API response cache")
    p_cache.add_argument("--clear", action="store_true", help="Delete all cached API responses")

    # --- serve ---
    p_serve = subparsers.add_parser(
        "serve",
        help="Start the interactive web application locally or in Hugging Face Spaces",
    )
    p_serve.add_argument("--host", default="127.0.0.1", help="Bind host (default: 127.0.0.1; use 0.0.0.0 in Docker)")
    p_serve.add_argument("--port", type=int, default=7860, help="Bind port (default: 7860)")
    p_serve.add_argument("--no-cache", action="store_true", help="Disable SQLite HTTP cache")
    p_serve.add_argument(
        "--verify-network",
        action="store_true",
        help="Verify outbound DNS resolution before binding port",
    )

    return parser


def _parse_batch_targets(raw_arg: str) -> list[str]:
    p = Path(raw_arg)
    if p.exists() and p.is_file():
        genes: list[str] = []
        for line in p.read_text(encoding="utf-8").splitlines():
            cleaned = line.split("#", 1)[0].strip()
            if cleaned:
                genes.extend(tok.strip() for tok in cleaned.split(",") if tok.strip())
        return genes
    return [tok.strip() for tok in raw_arg.split(",") if tok.strip()]


def main(argv: Optional[list[str]] = None) -> int:
    """CLI entrypoint."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "cache":
        http = CachedHttpClient()
        if args.clear:
            removed = http.clear_cache()
            print(f"Cleared {removed} cached API response(s) from {http.cache_path}")
        else:
            stats = http.cache_stats()
            print(f"Cache path    : {stats['path']}")
            print(f"Cached entries: {stats['entries']}")
        return 0

    if args.command == "check":
        http = CachedHttpClient(use_cache=not args.no_cache)
        runner = PreflightRunner(http=http, check_pams=not args.no_pam)
        try:
            report = runner.run(gene_query=args.gene, tissue=args.tissue, cell_line=args.cell_line)
        except Exception as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1

        if not args.quiet:
            sys.stdout.write(render_terminal(report, show_all_exons=args.all_exons))

        if args.html:
            out_html = Path(args.html)
            out_html.parent.mkdir(parents=True, exist_ok=True)
            out_html.write_text(render_html(report), encoding="utf-8")
            print(f"[+] Saved interactive HTML report -> {out_html}")

        if args.markdown:
            out_md = Path(args.markdown)
            out_md.parent.mkdir(parents=True, exist_ok=True)
            out_md.write_text(render_markdown(report), encoding="utf-8")
            print(f"[+] Saved Markdown dossier -> {out_md}")

        if args.json:
            out_json = Path(args.json)
            out_json.parent.mkdir(parents=True, exist_ok=True)
            out_json.write_text(json.dumps(report.to_dict(), indent=2) + "\n", encoding="utf-8")
            print(f"[+] Saved JSON report -> {out_json}")

        return 0

    if args.command == "batch":
        genes = _parse_batch_targets(args.targets)
        if not genes:
            print("Error: No target genes provided.", file=sys.stderr)
            return 1

        http = CachedHttpClient(use_cache=not args.no_cache)
        runner = PreflightRunner(http=http, check_pams=not args.no_pam)
        out_csv = Path(args.output)
        out_csv.parent.mkdir(parents=True, exist_ok=True)

        html_dir = Path(args.html_dir) if args.html_dir else None
        if html_dir:
            html_dir.mkdir(parents=True, exist_ok=True)

        rows_out: list[dict[str, object]] = []
        for idx, g in enumerate(genes, start=1):
            print(f"[{idx}/{len(genes)}] Auditing {g}...")
            try:
                rep = runner.run(gene_query=g, tissue=args.tissue, cell_line=args.cell_line)
            except Exception as exc:
                print(f"  [!] Failed on {g}: {exc}", file=sys.stderr)
                rows_out.append({"gene": g, "status": f"ERROR: {exc}"})
                continue

            coding = [e for e in rep.exons if e.is_coding]
            best_ex = max(coding, key=lambda e: e.priority_score) if coding else None
            top_par = rep.paralogs[0] if rep.paralogs else None
            rows_out.append(
                {
                    "gene": rep.gene_symbol,
                    "ensembl_id": rep.ensembl_id,
                    "status": "OK",
                    "risk_level": rep.modality.risk_level,
                    "recommended_modality": rep.modality.primary_modality,
                    "cell_line_chronos": rep.essentiality.matched_cell_line_chronos,
                    "pan_essential_fraction": rep.essentiality.depmap_pan_essential_fraction,
                    "gnomad_loeuf": rep.essentiality.gnomad_loeuf,
                    "tissue_tpm": rep.expression.requested_tissue_tpm,
                    "top_paralog": f"{top_par.gene_symbol} ({top_par.percent_identity:.0f}%, {top_par.buffering_risk})"
                    if top_par
                    else "",
                    "best_exon": f"Exon {best_ex.exon_number} ({best_ex.priority_score}/100)" if best_ex else "",
                    "best_crispick_region": best_ex.crispick_target_region if best_ex else "",
                    "chemical_probes": "; ".join(c.name for c in rep.validation.orthogonal_chemical_tools[:3]),
                }
            )
            if html_dir:
                (html_dir / f"{rep.gene_symbol}_preflight.html").write_text(render_html(rep), encoding="utf-8")

        fieldnames = [
            "gene",
            "ensembl_id",
            "status",
            "risk_level",
            "recommended_modality",
            "cell_line_chronos",
            "pan_essential_fraction",
            "gnomad_loeuf",
            "tissue_tpm",
            "top_paralog",
            "best_exon",
            "best_crispick_region",
            "chemical_probes",
        ]
        with out_csv.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for r in rows_out:
                writer.writerow(r)
        print(f"[+] Wrote batch summary ({len(rows_out)} target(s)) -> {out_csv}")
        return 0

    if args.command == "serve":
        try:
            run_web_server(
                host=args.host,
                port=args.port,
                use_cache=not args.no_cache,
                verify_network=args.verify_network,
            )
        except Exception as exc:
            print(f"Error starting server: {exc}", file=sys.stderr)
            return 1
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

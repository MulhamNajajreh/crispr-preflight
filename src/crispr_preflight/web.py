"""Zero-dependency local & Hugging Face Spaces web application server for crispr-preflight.

Runs on Python's standard-library ThreadingHTTPServer (default port 7860 to match
Hugging Face Spaces out of the box), serving an interactive web UI and JSON API.
"""

from __future__ import annotations

import json
import socket
import threading
import time
import urllib.parse
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Optional

from crispr_preflight import __version__
from crispr_preflight.cache import CachedHttpClient, NetworkFetchError
from crispr_preflight.engine.runner import PreflightRunner
from crispr_preflight.renderers.html_report import render_html
from crispr_preflight.renderers.markdown import render_markdown


WEB_APP_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>CRISPR-Preflight — Pre-Experiment Target &amp; Knockout Flight Check</title>
<style>
  :root {
    --bg: #f8fafc;
    --card: #ffffff;
    --fg: #0f172a;
    --muted: #64748b;
    --border: #e2e8f0;
    --primary: #0f766e;
    --primary-hover: #115e59;
    --primary-fg: #ffffff;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    background: var(--bg);
    color: var(--fg);
    line-height: 1.5;
  }
  .topbar {
    background: #0f172a;
    color: #f8fafc;
    padding: 16px 24px;
    border-bottom: 1px solid #1e293b;
  }
  .topbar-inner {
    max-width: 1180px;
    margin: 0 auto;
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between;
    gap: 14px;
  }
  .brand { font-size: 20px; font-weight: 800; letter-spacing: -0.02em; }
  .brand span { color: #2dd4bf; }
  .tagline { font-size: 12.5px; color: #94a3b8; }
  .mode-tabs {
    display: inline-flex;
    background: #1e293b;
    padding: 4px;
    border-radius: 10px;
    gap: 4px;
  }
  .mode-tab {
    border: none;
    background: transparent;
    color: #cbd5e1;
    font-size: 13px;
    font-weight: 700;
    padding: 7px 14px;
    border-radius: 7px;
    cursor: pointer;
  }
  .mode-tab:hover { color: #ffffff; }
  .mode-tab.active {
    background: #0f766e;
    color: #ffffff;
  }
  .wrap { max-width: 1180px; margin: 22px auto; padding: 0 20px 48px; }
  .card {
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 20px 24px;
    box-shadow: 0 2px 10px rgba(15, 23, 42, 0.04);
    margin-bottom: 18px;
  }
  .search-grid {
    display: grid;
    grid-template-columns: 1.08fr 1.22fr 0.95fr 1.05fr auto;
    gap: 12px;
    align-items: end;
  }
  .batch-grid {
    display: grid;
    grid-template-columns: 1.05fr 1.45fr 0.9fr 1fr auto;
    gap: 12px;
    align-items: end;
  }
  @media (max-width: 980px) {
    .search-grid, .batch-grid { grid-template-columns: 1fr; }
  }
  label {
    display: block;
    font-size: 12px;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--muted);
    margin-bottom: 6px;
  }
  input, select, textarea {
    width: 100%;
    padding: 9px 12px;
    font-size: 14px;
    font-family: inherit;
    border: 1px solid #cbd5e1;
    border-radius: 9px;
    background: #ffffff;
    color: var(--fg);
    transition: border-color 0.15s ease, box-shadow 0.15s ease;
  }
  input:not([type="checkbox"]), select {
    height: 40px;
  }
  input:hover, select:hover, textarea:hover {
    border-color: #94a3b8;
  }
  input:focus, select:focus, textarea:focus {
    outline: none;
    border-color: #0d9488;
    box-shadow: 0 0 0 3px rgba(20, 184, 166, 0.18);
  }
  select.nice-select {
    appearance: none;
    -webkit-appearance: none;
    padding-right: 34px;
    cursor: pointer;
    background-image: url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='16' height='16' viewBox='0 0 24 24' fill='none' stroke='%23475569' stroke-width='2.3' stroke-linecap='round' stroke-linejoin='round'%3E%3Cpolyline points='6 9 12 15 18 9'%3E%3C/polyline%3E%3C/svg%3E");
    background-repeat: no-repeat;
    background-position: right 10px center;
  }
  select:disabled {
    background-color: #f1f5f9;
    color: #94a3b8;
    cursor: not-allowed;
  }
  .combo-wrap {
    position: relative;
    width: 100%;
  }
  .combo-wrap input {
    padding-right: 58px;
  }
  .combo-actions {
    position: absolute;
    right: 6px;
    top: 50%;
    transform: translateY(-50%);
    display: inline-flex;
    align-items: center;
    gap: 2px;
  }
  .combo-icon-btn {
    border: none;
    background: transparent;
    color: #64748b;
    width: 24px;
    height: 26px;
    border-radius: 6px;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    font-size: 14px;
    padding: 0;
  }
  .combo-icon-btn:hover {
    background: #f1f5f9;
    color: #0f172a;
  }
  .combo-menu {
    position: absolute;
    top: calc(100% + 5px);
    left: 0;
    right: 0;
    min-width: 280px;
    max-height: 310px;
    overflow-y: auto;
    background: #ffffff;
    border: 1px solid #cbd5e1;
    border-radius: 11px;
    box-shadow: 0 12px 28px rgba(15, 23, 42, 0.14);
    z-index: 60;
    padding: 5px;
    display: none;
  }
  .combo-menu.open { display: block; }
  .combo-grp {
    font-size: 10.5px;
    font-weight: 800;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: #64748b;
    padding: 7px 9px 3px;
  }
  .combo-opt {
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 8px;
    width: 100%;
    text-align: left;
    border: none;
    background: transparent;
    padding: 7px 9px;
    border-radius: 7px;
    cursor: pointer;
    font-size: 13px;
    color: #0f172a;
  }
  .combo-opt:hover, .combo-opt.selected {
    background: #f0fdfa;
    color: #0f766e;
  }
  .combo-opt-sub {
    font-size: 11.5px;
    color: #64748b;
    font-weight: 500;
  }
  .btn-primary {
    padding: 10px 20px;
    height: 40px;
    font-size: 14px;
    font-weight: 700;
    background: var(--primary);
    color: var(--primary-fg);
    border: none;
    border-radius: 9px;
    cursor: pointer;
    white-space: nowrap;
  }
  .btn-primary:hover { background: var(--primary-hover); }
  .btn-primary:disabled { opacity: 0.6; cursor: wait; }
  .btn-secondary {
    padding: 6px 11px;
    font-size: 12px;
    font-weight: 600;
    background: #f1f5f9;
    color: #334155;
    border: 1px solid #cbd5e1;
    border-radius: 7px;
    cursor: pointer;
    white-space: nowrap;
  }
  .btn-secondary:hover { background: #e2e8f0; }
  .btn-accent {
    padding: 6px 12px;
    font-size: 12px;
    font-weight: 700;
    background: #f0fdfa;
    color: #0f766e;
    border: 1px solid #5eead4;
    border-radius: 7px;
    cursor: pointer;
    white-space: nowrap;
  }
  .btn-accent:hover { background: #ccfbf1; }
  .seg-bar {
    display: inline-flex;
    flex-wrap: wrap;
    background: #f1f5f9;
    border: 1px solid #cbd5e1;
    border-radius: 9px;
    padding: 3px;
    gap: 3px;
  }
  .seg-btn {
    border: none;
    background: transparent;
    color: #475569;
    font-size: 12px;
    font-weight: 700;
    padding: 5px 11px;
    border-radius: 6px;
    cursor: pointer;
  }
  .seg-btn.active {
    background: #ffffff;
    color: #0f172a;
    box-shadow: 0 1px 3px rgba(15, 23, 42, 0.1);
  }
  .chips { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 14px; align-items: center; }
  .chip-lbl { font-size: 12px; color: var(--muted); font-weight: 600; }
  .chip {
    font-size: 12px;
    padding: 4px 10px;
    border-radius: 999px;
    background: #f1f5f9;
    border: 1px solid #cbd5e1;
    color: #334155;
    cursor: pointer;
    font-weight: 600;
  }
  .chip:hover { background: #ccfbf1; border-color: #5eead4; color: #0f766e; }
  .status-box {
    padding: 15px 18px;
    border-radius: 10px;
    font-size: 14px;
    margin-bottom: 18px;
    display: none;
  }
  .status-loading { display: block; background: #f0fdfa; border: 1px solid #99f6e4; color: #0f766e; }
  .status-error { display: block; background: #fef2f2; border: 1px solid #fecaca; color: #b91c1c; }
  .verdict { border-radius: 12px; padding: 18px 22px; margin-bottom: 18px; border: 1px solid var(--border); }
  .verdict-HIGH { background: #fef2f2; border-color: #fecaca; color: #991b1b; }
  .verdict-MODERATE { background: #fffbeb; border-color: #fde68a; color: #92400e; }
  .verdict-LOW { background: #ecfdf5; border-color: #a7f3d0; color: #065f46; }
  .kpi-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; margin-bottom: 18px; }
  .kpi { background: var(--card); border: 1px solid var(--border); border-radius: 12px; padding: 14px 16px; }
  .kpi .lbl { font-size: 12px; color: var(--muted); font-weight: 600; text-transform: uppercase; }
  .kpi .val { font-size: 18px; font-weight: 800; margin-top: 4px; }
  .kpi .sub { font-size: 12px; color: var(--muted); margin-top: 2px; }
  h2 { font-size: 16px; margin: 0 0 12px 0; }
  .inspector {
    background: #f8fafc;
    border: 1px solid #cbd5e1;
    border-radius: 10px;
    padding: 14px 16px;
    margin-top: 12px;
  }
  .trap { border-radius: 9px; padding: 12px 16px; margin-bottom: 10px; background: #f8fafc; border: 1px solid var(--border); border-left-width: 4px; font-size: 13px; }
  .trap-HIGH { border-left-color: #ef4444; background: #fff5f5; }
  .trap-MODERATE { border-left-color: #f59e0b; background: #fffdf2; }
  .trap-INFO { border-left-color: #3b82f6; background: #f5f9ff; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  th, td { text-align: left; padding: 9px 10px; border-bottom: 1px solid var(--border); vertical-align: top; }
  th { font-size: 12px; color: var(--muted); background: #f8fafc; font-weight: 700; }
  tr.clickable-row { cursor: pointer; }
  tr.clickable-row:hover { background: #f1f5f9; }
  tr.selected-row { background: #f0fdfa; }
  .pill { display: inline-block; padding: 2px 8px; border-radius: 999px; font-size: 11.5px; font-weight: 700; }
  .pill-RECOMMENDED { background: #d1fae5; color: #065f46; }
  .pill-CAUTION { background: #fef3c7; color: #92400e; }
  .pill-AVOID { background: #fee2e2; color: #991b1b; }
  .pill-INFO { background: #e0f2fe; color: #075985; }
  .dom-legend { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 10px; padding: 10px 12px; background: #eef2ff; border: 1px solid #c7d2fe; border-radius: 9px; align-items: center; }
  .dom-chip { display: inline-flex; align-items: center; gap: 6px; background: #ffffff; border: 1px solid #a5b4fc; border-radius: 7px; padding: 4px 9px; font-size: 12px; color: #1e1b4b; cursor: pointer; }
  .dom-chip:hover { background: #e0e7ff; border-color: #6366f1; }
  .dom-badge { background: #4f46e5; color: #ffffff; font-weight: 800; font-size: 11px; padding: 1px 6px; border-radius: 4px; }
  .pam-badge { background: #dbeafe; color: #1e40af; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 700; font-size: 12px; padding: 2px 5px; border-radius: 4px; }
  .offtarget-note { background: #eff6ff; border: 1px solid #bfdbfe; color: #1e3a8a; border-radius: 9px; padding: 10px 14px; font-size: 12.5px; margin-bottom: 12px; display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 10px; }
  .two-col { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  @media (max-width: 820px) { .two-col { grid-template-columns: 1fr; } }
  code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 12px; background: #f1f5f9; padding: 2px 5px; border-radius: 4px; }
  .ccle-tag { display: inline-block; margin-top: 5px; padding: 2px 7px; border-radius: 6px; font-size: 11.5px; font-weight: 600; }
  .ccle-alert { background: #fee2e2; color: #991b1b; border: 1px solid #fecaca; }
  .ccle-warn { background: #fef3c7; color: #92400e; border: 1px solid #fde68a; }
  .ccle-ok { background: #f1f5f9; color: #334155; border: 1px solid #cbd5e1; }
</style>
</head>
<body>
  <div class="topbar">
    <div class="topbar-inner">
      <div>
        <div class="brand">crispr-<span>preflight</span></div>
        <div class="tagline">Pre-experiment biological flight check across Human, Mouse, Rat, Zebrafish, Drosophila &amp; C. elegans</div>
      </div>
      <div class="mode-tabs" role="tablist">
        <button type="button" id="tab-single" class="mode-tab active" onclick="switchAppMode('SINGLE')">Single-Gene Flight Check</button>
        <button type="button" id="tab-batch" class="mode-tab" onclick="switchAppMode('BATCH')">Batch Hit-List Triage (5–30 Genes)</button>
      </div>
    </div>
  </div>

  <div class="wrap">
    <!-- ==================== MODE 1: SINGLE-GENE FLIGHT CHECK ==================== -->
    <div id="single-mode-view">
      <div class="card">
        <form id="query-form" onsubmit="runCheck(event)">
          <div class="search-grid">
            <div>
              <label for="species-input">Organism / Genome</label>
              <select id="species-input" class="nice-select" onchange="onSpeciesChange('single')">
                <option value="human" selected>Human — H. sapiens (GRCh38)</option>
                <option value="mouse">Mouse — M. musculus (GRCm39)</option>
                <option value="rat">Rat — R. norvegicus (GRCr8)</option>
                <option value="zebrafish">Zebrafish — D. rerio (GRCz11)</option>
                <option value="fly">Fruit Fly — D. melanogaster (BDGP6.46)</option>
                <option value="worm">Roundworm — C. elegans (WBcel235)</option>
              </select>
            </div>
            <div>
              <label for="gene-input" id="gene-input-lbl">Gene Symbol, Alias, or ID</label>
              <input id="gene-input" type="text" placeholder="Type a gene (e.g. STAT3, TP53, BRD4, KRAS...)" value="" autofocus required />
            </div>
            <div>
              <label for="tissue-input" id="tissue-input-lbl">Target Tissue / Lineage (Optional)</label>
              <select id="tissue-input" class="nice-select">
                <option value="" id="single-tissue-opt0" selected>Any / Pan-Tissue (Peak GTEx)</option>
                <optgroup label="Blood &amp; Immune">
                  <option value="blood">Blood / Myeloid / Lymphoid</option>
                </optgroup>
                <optgroup label="Epithelial &amp; Solid Organs">
                  <option value="lung">Lung</option>
                  <option value="colon">Colon / Gut</option>
                  <option value="breast">Breast</option>
                  <option value="liver">Liver</option>
                  <option value="pancreas">Pancreas</option>
                  <option value="kidney">Kidney</option>
                  <option value="prostate">Prostate</option>
                  <option value="ovary">Ovary</option>
                  <option value="heart">Heart</option>
                </optgroup>
                <optgroup label="Neural, Skin &amp; Stromal">
                  <option value="brain">Brain / Neural</option>
                  <option value="skin">Skin / Melanocyte</option>
                  <option value="fibroblast">Cultured Fibroblasts</option>
                </optgroup>
              </select>
            </div>
            <div>
              <label for="cell-input" id="cell-input-lbl">Cell Line / Model (Optional)</label>
              <div class="combo-wrap" id="single-combo-wrap">
                <input id="cell-input" type="text" autocomplete="off" placeholder="Optional — browse or type model..." value=""
                  onfocus="openCellMenu('single', true)" oninput="onCellInput('single')" onkeydown="onCellKeyDown(event, 'single')" />
                <div class="combo-actions">
                  <button type="button" id="single-cell-clear" class="combo-icon-btn" style="display:none;" onclick="clearCellInput('single')" title="Clear cell line">&times;</button>
                  <button type="button" class="combo-icon-btn" onclick="toggleCellMenu(event, 'single')" title="Browse common models">&#9662;</button>
                </div>
                <div id="single-cell-menu" class="combo-menu"></div>
              </div>
            </div>
            <div>
              <button id="run-btn" type="submit" class="btn-primary">Run Flight Check</button>
            </div>
          </div>
        </form>

        <div class="chips" id="example-chips-bar"></div>
      </div>

      <div id="status-box" class="status-box"></div>
    </div>

    <!-- ==================== MODE 2: BATCH HIT-LIST TRIAGE ==================== -->
    <div id="batch-mode-view" style="display:none;">
      <div class="card">
        <div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:10px;margin-bottom:12px;">
          <div>
            <h2 style="margin:0;font-size:18px;">Batch Hit-List Triage (Paste 5–30 Screen or RNA-seq Candidate Genes)</h2>
            <div style="font-size:13px;color:var(--muted);">Compares modality risk, essentiality, cell-line genetics, paralog traps, and top safe exons side-by-side — with full per-gene reports right below.</div>
          </div>
          <button type="button" class="btn-secondary" id="demo-batch-btn" onclick="loadDemoBatch()">Load Demo Screen Hits</button>
        </div>

        <div class="batch-grid">
          <div>
            <label for="batch-species">Organism / Genome</label>
            <select id="batch-species" class="nice-select" onchange="onSpeciesChange('batch')">
              <option value="human" selected>Human — H. sapiens (GRCh38)</option>
              <option value="mouse">Mouse — M. musculus (GRCm39)</option>
              <option value="rat">Rat — R. norvegicus (GRCr8)</option>
              <option value="zebrafish">Zebrafish — D. rerio (GRCz11)</option>
              <option value="fly">Fruit Fly — D. melanogaster (BDGP6.46)</option>
              <option value="worm">Roundworm — C. elegans (WBcel235)</option>
            </select>
          </div>
          <div>
            <label for="batch-genes">Gene Symbols (comma, space, or newline separated — up to 30)</label>
            <textarea id="batch-genes" rows="2" placeholder="Paste gene list here (e.g. BRD4, SMARCA4, PLK1, SOX2, STAT3, KRAS)..."></textarea>
          </div>
          <div>
            <label for="batch-tissue">Target Tissue (Optional)</label>
            <select id="batch-tissue" class="nice-select">
              <option value="" id="batch-tissue-opt0" selected>Any / Pan-Tissue (Peak GTEx)</option>
              <optgroup label="Blood &amp; Immune">
                <option value="blood">Blood / Myeloid / Lymphoid</option>
              </optgroup>
              <optgroup label="Epithelial &amp; Solid Organs">
                <option value="lung">Lung</option>
                <option value="colon">Colon / Gut</option>
                <option value="breast">Breast</option>
                <option value="liver">Liver</option>
                <option value="pancreas">Pancreas</option>
                <option value="kidney">Kidney</option>
                <option value="prostate">Prostate</option>
                <option value="ovary">Ovary</option>
                <option value="heart">Heart</option>
              </optgroup>
              <optgroup label="Neural, Skin &amp; Stromal">
                <option value="brain">Brain / Neural</option>
                <option value="skin">Skin / Melanocyte</option>
                <option value="fibroblast">Cultured Fibroblasts</option>
              </optgroup>
            </select>
          </div>
          <div>
            <label for="batch-cell">Cell Line / Model (Optional)</label>
            <div class="combo-wrap" id="batch-combo-wrap">
              <input id="batch-cell" type="text" autocomplete="off" value="" placeholder="Optional — browse or type model..."
                onfocus="openCellMenu('batch', true)" oninput="onCellInput('batch')" onkeydown="onCellKeyDown(event, 'batch')" />
              <div class="combo-actions">
                <button type="button" id="batch-cell-clear" class="combo-icon-btn" style="display:none;" onclick="clearCellInput('batch')" title="Clear cell line">&times;</button>
                <button type="button" class="combo-icon-btn" onclick="toggleCellMenu(event, 'batch')" title="Browse common models">&#9662;</button>
              </div>
              <div id="batch-cell-menu" class="combo-menu"></div>
            </div>
          </div>
          <div>
            <button type="button" id="batch-run-btn" class="btn-primary" onclick="runBatchTriage()">Triage Gene List</button>
          </div>
        </div>
        <div id="batch-progress" style="font-size:13px;color:var(--primary);font-weight:600;margin-top:12px;display:none;"></div>
      </div>

      <div id="batch-results-card" class="card" style="display:none;">
        <div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:8px;margin-bottom:12px;">
          <div>
            <h2 style="margin:0;">Hit-List Comparison Matrix</h2>
            <div style="font-size:12px;color:var(--muted);">Click any gene row to inspect its full interactive exon map, sgRNAs &amp; PCR primers right below on this page.</div>
          </div>
          <button type="button" class="btn-accent" onclick="exportBatchCsv()">Export Triage Summary (.csv)</button>
        </div>
        <div style="overflow-x:auto;">
          <table>
            <thead>
              <tr>
                <th>Gene</th>
                <th>Risk</th>
                <th>Recommended Modality</th>
                <th>Cell Fitness (Chronos)</th>
                <th>CCLE Cell-Line Genetics</th>
                <th>Tissue TPM</th>
                <th>Top Paralog</th>
                <th>Best Exon &amp; Top Spacer</th>
                <th>Full Report</th>
              </tr>
            </thead>
            <tbody id="batch-tbody"></tbody>
          </table>
        </div>
      </div>
    </div>

    <!-- ==================== SHARED FULL INTERACTIVE GENE REPORT ==================== -->
    <div id="results-area" style="display:none;">
      <div id="batch-report-nav" class="card" style="display:none;background:#f0fdfa;border-color:#99f6e4;padding:12px 18px;"></div>

      <!-- Gene Header & Export / Order Bar -->
      <div class="card" style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:12px;">
        <div>
          <div style="font-size:22px;font-weight:800;" id="res-title"></div>
          <div style="font-size:13px;color:var(--muted);" id="res-sub"></div>
        </div>
        <div style="display:flex;flex-wrap:wrap;gap:7px;align-items:center;">
          <button type="button" class="btn-accent" onclick="exportIdtOrderCsv()">Export IDT Oligo Order (.csv)</button>
          <button type="button" id="basket-btn" class="btn-accent" style="display:none;" onclick="exportBasketOnlyCsv()"></button>
          <button type="button" id="basket-clear-btn" class="btn-secondary" style="display:none;" onclick="clearOrderBasket()" title="Clear pinned guides">Clear Basket</button>
          <button type="button" class="btn-secondary" onclick="copyShareableLink(this)">Copy Share Link</button>
          <button type="button" class="btn-secondary" onclick="downloadReport('html')">HTML (Benchling)</button>
          <button type="button" class="btn-secondary" onclick="downloadReport('md')">Markdown</button>
          <button type="button" class="btn-secondary" onclick="downloadReport('json')">JSON</button>
        </div>
      </div>

      <!-- Modality Verdict -->
      <div id="res-verdict" class="verdict"></div>

      <!-- KPIs -->
      <div class="kpi-grid">
        <div class="kpi">
          <div class="lbl" id="lbl-kpi-1">DepMap Fitness &amp; CCLE Genetics</div>
          <div class="val" id="k-chronos"></div>
          <div class="sub" id="k-chronos-sub"></div>
          <div id="k-ccle-sub"></div>
        </div>
        <div class="kpi"><div class="lbl" id="lbl-kpi-2">Human LoF Constraint (gnomAD)</div><div class="val" id="k-loeuf"></div><div class="sub" id="k-pli-sub"></div></div>
        <div class="kpi"><div class="lbl" id="lbl-kpi-3">GTEx Baseline Expression</div><div class="val" id="k-tpm"></div><div class="sub" id="k-loc"></div></div>
        <div class="kpi"><div class="lbl">Audited Coding Isoforms</div><div class="val" id="k-iso"></div><div class="sub" id="k-exons"></div></div>
      </div>

      <!-- SVG Architecture -->
      <div class="card">
        <div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:8px;margin-bottom:8px;">
          <div>
            <h2 style="margin:0;">1. Interactive Exon, Functional Domain &amp; Epitope Map (5&prime; &rarr; 3&prime;)</h2>
            <div style="font-size:12px;color:var(--muted);">Click any colored exon block in Track 1 to inspect its score math and exon-specific guides.</div>
          </div>
          <div style="font-size:12px;display:flex;gap:12px;">
            <span><b style="color:#059669;">&#9632;</b> Recommended (&ge;75)</span>
            <span><b style="color:#d97706;">&#9632;</b> Caution (45–74)</span>
            <span><b style="color:#e11d48;">&#9632;</b> Avoid (&lt;45)</span>
          </div>
        </div>
        <div id="svg-wrap" style="overflow-x:auto;"></div>
        <div id="domain-legend"></div>
        <div id="exon-inspector" class="inspector"></div>
      </div>

      <!-- Trap Warnings -->
      <div class="card">
        <h2>2. Biological &amp; Cell-Line Failure-Mode Alerts</h2>
        <div id="traps-wrap"></div>
      </div>

      <!-- Exon Table -->
      <div class="card">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:10px;">
          <h2 style="margin:0;">3. Coding Exon Targetability Audit</h2>
          <button type="button" class="btn-secondary" id="filter-btn" onclick="toggleExonFilter()">Showing: All Coding Exons</button>
        </div>
        <div style="overflow-x:auto;">
          <table>
            <thead>
              <tr>
                <th>Exon</th><th>Score</th><th>Tier</th><th>AA Span</th><th>Length</th><th>Isoforms</th><th>NGG Sites</th><th>NMD Status</th><th id="exon-coord-th">Target Region (GRCh38)</th>
              </tr>
            </thead>
            <tbody id="exon-tbody"></tbody>
          </table>
        </div>
      </div>

      <!-- Candidate SpCas9 KO & CRISPRi Guides + SNP Shield + TIDE/ICE Primers -->
      <div class="card">
        <div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:10px;margin-bottom:10px;">
          <div>
            <h2 style="margin:0;">4. SpCas9 sgRNA Spacers, Strain/SNP Shield &amp; TIDE/ICE Genotyping Primers</h2>
            <div style="font-size:12px;color:var(--muted);" id="guides-sub">Includes U6 Golden Gate oligos (<code>BsmBI</code>/<code>BbsI</code>), paired 430–740 bp TIDE/ICE PCR primers, and background gnomAD common SNP screening (MAF &ge;1%).</div>
          </div>
          <div style="display:flex;flex-wrap:wrap;gap:8px;align-items:center;">
            <div class="seg-bar">
              <button type="button" id="scope-panel" class="seg-btn active" onclick="setGuideScope('PANEL')">Top KO Panel</button>
              <button type="button" id="scope-exon" class="seg-btn" onclick="setGuideScope('EXON')">Selected Exon KO</button>
              <button type="button" id="scope-crispri" class="seg-btn" onclick="setGuideScope('CRISPRI')">CRISPRi Promoter (TSS)</button>
            </div>
            <button type="button" class="btn-secondary" onclick="copyGuidesFasta(this)">Copy FASTA</button>
          </div>
        </div>
        <div class="offtarget-note">
          <div style="max-width:740px;" id="offtarget-assembly-note">
            <strong>Before ordering oligos — 10-second genome-wide off-target check:</strong> Local rules rank on-target quality &amp; biological safety; paste your chosen 2–3 spacers into <strong>Cas-OFFinder</strong> or <strong>Broad CRISPick</strong> to verify zero 0–2 mismatch off-target sites across GRCh38.
          </div>
          <div style="display:flex;gap:8px;white-space:nowrap;">
            <a href="http://www.rgenome.net/cas-offinder/" target="_blank" rel="noopener" class="btn-secondary" style="text-decoration:none;">Cas-OFFinder &nearr;</a>
            <a href="https://portals.broadinstitute.org/gppx/crispick/public" target="_blank" rel="noopener" class="btn-secondary" style="text-decoration:none;">Broad CRISPick &nearr;</a>
          </div>
        </div>
        <div style="overflow-x:auto;" id="guides-wrap"></div>
      </div>

      <!-- Paralogs & Validation -->
      <div class="two-col">
        <div class="card">
          <div style="display:flex;justify-content:space-between;align-items:baseline;gap:8px;margin-bottom:8px;">
            <h2 style="margin:0;" id="paralog-hdr">5. Paralog Buffering Check</h2>
            <span style="font-size:11.5px;color:var(--muted);">Click a paralog to pin partner guides for Dual-KO</span>
          </div>
          <div id="paralogs-wrap"></div>
        </div>
        <div class="card">
          <h2>6. Wet-Lab Validation &amp; Orthogonal Chemical Tools</h2>
          <div id="val-wrap"></div>
        </div>
      </div>
    </div>

    <footer style="margin-top:28px;padding-top:14px;border-top:1px solid var(--border);font-size:12px;color:var(--muted);text-align:center;line-height:1.6;">
      <div><strong>crispr-preflight</strong> &bull; Developed by <strong>Mulham Najajreh</strong> &bull; Open-source under the MIT License</div>
      <div>This work was supported in part by complimentary access to Google Antigravity and Gemini 4 Argon to accelerate scientific workflows.</div>
    </footer>
  </div>

<script>
let APP_MODE = 'SINGLE'; // 'SINGLE' or 'BATCH'
let CURRENT = null;
let SINGLE_PAYLOAD = null;
let SELECTED_EXON = null;
let ONLY_REC = false;
let GUIDE_SCOPE = 'PANEL'; // 'PANEL', 'EXON', or 'CRISPRI'
let ACTIVE_REQ_ID = 0;
let ACTIVE_SNP_REQ = 0;
let SNP_CACHE = {}; // guide_id -> { status, summary, variants }
let ORDER_BASKET = []; // pinned guides across genes (e.g. target + paralog)
let BATCH_ROWS = []; // [{ gene, payload, error }]
let SELECTED_BATCH_IDX = -1;

const SPECIES_META = {
  human: {
    label: 'Human (GRCh38)',
    assembly: 'GRCh38',
    placeholder: 'Type a gene (e.g. STAT3, TP53, BRD4, KRAS...)',
    batchPlaceholder: 'Paste gene list here (e.g. BRD4, SMARCA4, PLK1, SOX2, STAT3, KRAS)...',
    demoBatch: 'BRD4, SMARCA4, PLK1, SOX2, STAT3, KRAS',
    cellSub: 'Pan-cancer median',
    customHint: 'Custom DepMap / CCLE Cell Line (1,739 supported)'
  },
  mouse: {
    label: 'Mouse (GRCm39)',
    assembly: 'GRCm39',
    placeholder: 'Type a mouse gene (e.g. Trp53, Smarca4, Plk1, Sox2...)',
    batchPlaceholder: 'Paste mouse genes (e.g. Trp53, Smarca4, Plk1, Sox2, Brd4, Stat3)...',
    demoBatch: 'Trp53, Smarca4, Plk1, Sox2, Brd4, Stat3',
    cellSub: 'General mouse model',
    customHint: 'Custom Mouse Cell Line / Strain'
  },
  rat: {
    label: 'Rat (GRCr8)',
    assembly: 'GRCr8',
    placeholder: 'Type a rat gene (e.g. Tp53, Smarca4, Stat3, Kras...)',
    batchPlaceholder: 'Paste rat genes (e.g. Tp53, Smarca4, Stat3, Kras)...',
    demoBatch: 'Tp53, Smarca4, Stat3, Kras',
    cellSub: 'General rat model',
    customHint: 'Custom Rat Cell Line / Strain'
  },
  zebrafish: {
    label: 'Zebrafish (GRCz11)',
    assembly: 'GRCz11',
    placeholder: 'Type a zebrafish gene (e.g. tp53, smarca4a, sox2, stat3...)',
    batchPlaceholder: 'Paste zebrafish genes (e.g. tp53, smarca4a, sox2, stat3)...',
    demoBatch: 'tp53, smarca4a, sox2, stat3',
    cellSub: 'F0 embryo / general model',
    customHint: 'Custom Zebrafish Line / Stage'
  },
  fly: {
    label: 'Drosophila (BDGP6.46)',
    assembly: 'BDGP6.46',
    placeholder: 'Type a fly gene (e.g. Stat92E, brm, Ras85D, N...)',
    batchPlaceholder: 'Paste fly genes (e.g. Stat92E, brm, Ras85D)...',
    demoBatch: 'Stat92E, brm, Ras85D',
    cellSub: 'S2R+ / germline stock',
    customHint: 'Custom Drosophila Cell Line / Stock'
  },
  worm: {
    label: 'C. elegans (WBcel235)',
    assembly: 'WBcel235',
    placeholder: 'Type a C. elegans gene (e.g. cep-1, swsn-4, let-60...)',
    batchPlaceholder: 'Paste C. elegans genes (e.g. cep-1, swsn-4, let-60)...',
    demoBatch: 'cep-1, swsn-4, let-60',
    cellSub: 'N2 Bristol / RNP injection',
    customHint: 'Custom C. elegans Strain'
  }
};

const EXAMPLES_BY_SPECIES = {
  human: [
    { gene: 'BRD4', tissue: 'lung', cell: 'A549', label: 'BRD4 (Isoforms + Probes)' },
    { gene: 'SMARCA4', tissue: 'lung', cell: 'A549', label: 'SMARCA4 (Already LoF in A549 + SMARCA2 Trap!)' },
    { gene: 'PLK1', tissue: 'lung', cell: 'A549', label: 'PLK1 (Pan-Essential \u2192 CRISPRi)' },
    { gene: 'SOX2', tissue: 'lung', cell: 'A549', label: 'SOX2 (1-Exon Intronless NMD Escape)' },
    { gene: 'STAT3', tissue: 'blood', cell: 'K562', label: 'STAT3 (Blood / K562)' },
    { gene: 'KRAS', tissue: 'colon', cell: 'HCT116', label: 'KRAS (Endogenous Driver in HCT116)' }
  ],
  mouse: [
    { gene: 'Trp53', tissue: '', cell: 'NIH/3T3', label: 'Trp53 (Mouse p53 \u2022 GRCm39)' },
    { gene: 'Smarca4', tissue: '', cell: 'B16-F10', label: 'Smarca4 (Smarca2 Paralog Trap!)' },
    { gene: 'Plk1', tissue: '', cell: 'mESC (E14)', label: 'Plk1 (Essential \u2192 CRISPRi)' },
    { gene: 'Sox2', tissue: '', cell: 'mESC (E14)', label: 'Sox2 (1-Exon Intronless)' },
    { gene: 'Brd4', tissue: '', cell: 'RAW 264.7', label: 'Brd4 (Isoforms + BET Probes)' },
    { gene: 'Stat3', tissue: '', cell: 'RAW 264.7', label: 'Stat3 (RAW 264.7 Myeloid)' }
  ],
  rat: [
    { gene: 'Tp53', tissue: '', cell: 'PC-12', label: 'Tp53 (Rat p53 \u2022 GRCr8)' },
    { gene: 'Smarca4', tissue: '', cell: 'H9c2', label: 'Smarca4 (Smarca2 Paralog Trap)' },
    { gene: 'Stat3', tissue: '', cell: 'NRK', label: 'Stat3 (Rat STAT3)' },
    { gene: 'Kras', tissue: '', cell: 'PC-12', label: 'Kras (Rat KRAS)' }
  ],
  zebrafish: [
    { gene: 'tp53', tissue: '', cell: 'AB / Tubingen Embryo (F0)', label: 'tp53 (Zebrafish GRCz11)' },
    { gene: 'smarca4a', tissue: '', cell: 'AB / Tubingen Embryo (F0)', label: 'smarca4a (Teleost Paralog Trap!)' },
    { gene: 'sox2', tissue: '', cell: 'AB / Tubingen Embryo (F0)', label: 'sox2 (1-Exon Intronless)' },
    { gene: 'stat3', tissue: '', cell: 'AB / Tubingen Embryo (F0)', label: 'stat3 (Zebrafish stat3)' }
  ],
  fly: [
    { gene: 'Stat92E', tissue: '', cell: 'S2R+', label: 'Stat92E (JAK/STAT \u2022 BDGP6.46)' },
    { gene: 'brm', tissue: '', cell: 'S2R+', label: 'brm (Brahma SWI/SNF)' },
    { gene: 'Ras85D', tissue: '', cell: 'Kc167', label: 'Ras85D (Essential \u2192 CRISPRi/RNAi)' }
  ],
  worm: [
    { gene: 'cep-1', tissue: '', cell: 'N2 Bristol Gonad (RNP)', label: 'cep-1 (Worm p53 \u2022 WBcel235)' },
    { gene: 'swsn-4', tissue: '', cell: 'N2 Bristol Gonad (RNP)', label: 'swsn-4 (SWI/SNF ATPase)' },
    { gene: 'let-60', tissue: '', cell: 'N2 Bristol Gonad (RNP)', label: 'let-60 (Ras Ortholog)' }
  ]
};

const CELL_CATALOG_BY_SPECIES = {
  human: [
    { group: 'Common Lab Workhorses', name: 'HEK293T', sub: 'Embryonic Kidney' },
    { group: 'Common Lab Workhorses', name: 'HeLa', sub: 'Cervical Carcinoma' },
    { group: 'Common Lab Workhorses', name: 'U2OS', sub: 'Osteosarcoma' },
    { group: 'Common Lab Workhorses', name: 'RPE1', sub: 'Retinal Epithelial (Diploid)' },
    { group: 'Common Lab Workhorses', name: 'IMR-90', sub: 'Lung Fibroblast' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'K562', sub: 'CML Myeloid' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'Jurkat', sub: 'T-Cell Leukemia (T-ALL)' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'THP-1', sub: 'AML Monocytic' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'MOLM-13', sub: 'AML (FLT3-ITD)' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'MV4-11', sub: 'AML Macrophage' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'NALM-6', sub: 'B-Cell Precursor ALL' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'Raji', sub: 'Burkitt B-Cell Lymphoma' },
    { group: 'Blood / Leukemia / Lymphoma', name: 'HL-60', sub: 'Promyelocytic Leukemia' },
    { group: 'Lung', name: 'A549', sub: 'NSCLC Adenocarcinoma' },
    { group: 'Lung', name: 'NCI-H1299', sub: 'NSCLC (p53-null)' },
    { group: 'Lung', name: 'NCI-H1975', sub: 'NSCLC (EGFR L858R/T790M)' },
    { group: 'Lung', name: 'PC-9', sub: 'NSCLC (EGFR exon19del)' },
    { group: 'Lung', name: 'HCC827', sub: 'NSCLC Adenocarcinoma' },
    { group: 'Lung', name: 'Calu-3', sub: 'Lung Airway Epithelial' },
    { group: 'Colon & GI', name: 'HCT116', sub: 'Colorectal (MSI-H, KRAS G13D)' },
    { group: 'Colon & GI', name: 'HT-29', sub: 'Colorectal (BRAF V600E)' },
    { group: 'Colon & GI', name: 'SW480', sub: 'Colorectal Adenocarcinoma' },
    { group: 'Colon & GI', name: 'DLD-1', sub: 'Colorectal Adenocarcinoma' },
    { group: 'Colon & GI', name: 'Caco-2', sub: 'Intestinal Epithelial Barrier' },
    { group: 'Colon & GI', name: 'LoVo', sub: 'Colorectal Adenocarcinoma' },
    { group: 'Breast', name: 'MCF7', sub: 'Luminal ER+ Breast' },
    { group: 'Breast', name: 'MDA-MB-231', sub: 'Triple-Negative Breast (TNBC)' },
    { group: 'Breast', name: 'SK-BR-3', sub: 'HER2-Amplified Breast' },
    { group: 'Breast', name: 'T-47D', sub: 'Luminal ER+ Breast' },
    { group: 'Breast', name: 'BT-474', sub: 'Luminal HER2+ Breast' },
    { group: 'Liver & Pancreas', name: 'HepG2', sub: 'Hepatoblastoma / Liver' },
    { group: 'Liver & Pancreas', name: 'Huh-7', sub: 'Hepatocellular Carcinoma' },
    { group: 'Liver & Pancreas', name: 'PANC-1', sub: 'Pancreatic Ductal (PDAC)' },
    { group: 'Liver & Pancreas', name: 'MIA PaCa-2', sub: 'Pancreatic Ductal (PDAC)' },
    { group: 'Liver & Pancreas', name: 'AsPC-1', sub: 'Pancreatic Adenocarcinoma' },
    { group: 'Brain, Skin, Prostate & Ovary', name: 'U-87 MG', sub: 'Glioblastoma' },
    { group: 'Brain, Skin, Prostate & Ovary', name: 'SH-SY5Y', sub: 'Neuroblastoma' },
    { group: 'Brain, Skin, Prostate & Ovary', name: 'A375', sub: 'Melanoma (BRAF V600E)' },
    { group: 'Brain, Skin, Prostate & Ovary', name: 'PC-3', sub: 'Prostate Carcinoma' },
    { group: 'Brain, Skin, Prostate & Ovary', name: 'LNCaP', sub: 'Androgen-Sensitive Prostate' },
    { group: 'Brain, Skin, Prostate & Ovary', name: 'OVCAR-8', sub: 'High-Grade Serous Ovarian' },
    { group: 'Brain, Skin, Prostate & Ovary', name: '786-O', sub: 'Clear-Cell Renal (VHL-null)' }
  ],
  mouse: [
    { group: 'Fibroblast & Stem Cell Workhorses', name: 'NIH/3T3', sub: 'Embryonic Fibroblast' },
    { group: 'Fibroblast & Stem Cell Workhorses', name: 'MEF', sub: 'Mouse Embryonic Fibroblast' },
    { group: 'Fibroblast & Stem Cell Workhorses', name: 'mESC (E14)', sub: 'Pluripotent Embryonic Stem Cell' },
    { group: 'Immune & Hematopoietic', name: 'RAW 264.7', sub: 'Macrophage / Myeloid' },
    { group: 'Immune & Hematopoietic', name: 'BV-2', sub: 'Microglia' },
    { group: 'Immune & Hematopoietic', name: 'EL4', sub: 'T-Cell Lymphoma' },
    { group: 'Immune & Hematopoietic', name: 'Ba/F3', sub: 'IL-3 Dependent Pro-B' },
    { group: 'Syngeneic Oncology Models', name: 'B16-F10', sub: 'Melanoma (C57BL/6)' },
    { group: 'Syngeneic Oncology Models', name: 'MC38', sub: 'Colon Adenocarcinoma (C57BL/6)' },
    { group: 'Syngeneic Oncology Models', name: 'CT26', sub: 'Colon Carcinoma (BALB/c)' },
    { group: 'Syngeneic Oncology Models', name: '4T1', sub: 'Mammary Carcinoma (BALB/c)' },
    { group: 'Muscle & Neural', name: 'C2C12', sub: 'Myoblast' },
    { group: 'Muscle & Neural', name: 'Neuro-2a', sub: 'Neuroblastoma' }
  ],
  rat: [
    { group: 'Rat Cell Models', name: 'PC-12', sub: 'Adrenal Pheochromocytoma / Neuronal' },
    { group: 'Rat Cell Models', name: 'H9c2', sub: 'Embryonic Cardiomyoblast' },
    { group: 'Rat Cell Models', name: 'NRK', sub: 'Normal Rat Kidney Epithelial/Fibroblast' },
    { group: 'Rat Cell Models', name: 'INS-1', sub: 'Pancreatic Beta-Cell' },
    { group: 'Rat Cell Models', name: 'C6', sub: 'Glioma' }
  ],
  zebrafish: [
    { group: 'Zebrafish In Vivo & Cell Lines', name: 'AB / Tubingen Embryo (F0)', sub: '1-Cell Zygote Cas9 RNP Injection' },
    { group: 'Zebrafish In Vivo & Cell Lines', name: 'TL Embryo', sub: '1-Cell Zygote Injection' },
    { group: 'Zebrafish In Vivo & Cell Lines', name: 'ZF4', sub: 'Embryonic Fibroblast Line' },
    { group: 'Zebrafish In Vivo & Cell Lines', name: 'PAC2', sub: 'Embryonic Fibroblast Line' }
  ],
  fly: [
    { group: 'Drosophila Tissue Culture & Germline', name: 'S2R+', sub: 'Adherent Schneider Embryonic' },
    { group: 'Drosophila Tissue Culture & Germline', name: 'Kc167', sub: 'Embryonic Hemocyte-like' },
    { group: 'Drosophila Tissue Culture & Germline', name: 'nos-Cas9 Germline', sub: 'Transgenic Germline Injection' }
  ],
  worm: [
    { group: 'C. elegans In Vivo Models', name: 'N2 Bristol Gonad (RNP)', sub: 'Syncytial Gonad Cas9 RNP + dpy-10' }
  ]
};

function getScopeSpecies(scope) {
  const el = document.getElementById(scope === 'batch' ? 'batch-species' : 'species-input');
  const val = (el && el.value ? el.value.trim().toLowerCase() : 'human');
  return SPECIES_META[val] ? val : 'human';
}

function getCatalogForScope(scope) {
  return CELL_CATALOG_BY_SPECIES[getScopeSpecies(scope)] || CELL_CATALOG_BY_SPECIES.human;
}

function renderExampleChips() {
  const sp = getScopeSpecies('single');
  const list = EXAMPLES_BY_SPECIES[sp] || EXAMPLES_BY_SPECIES.human;
  const bar = document.getElementById('example-chips-bar');
  if (!bar) return;
  bar.innerHTML = `<span class="chip-lbl">Try classic ${esc(SPECIES_META[sp].label.split(' ')[0])} examples:</span>` +
    list.map((ex, i) => `<button type="button" class="chip" onclick="runExampleChip(${i})">${esc(ex.label)}</button>`).join('');
}

function runExampleChip(idx) {
  const sp = getScopeSpecies('single');
  const list = EXAMPLES_BY_SPECIES[sp] || EXAMPLES_BY_SPECIES.human;
  const item = list[idx];
  if (!item) return;
  quickRun(item.gene, item.tissue, item.cell, sp);
}

function onSpeciesChange(scope, clearInputs = true) {
  const sp = getScopeSpecies(scope);
  const meta = SPECIES_META[sp];
  const isHuman = sp === 'human';
  if (scope === 'batch') {
    const tSel = document.getElementById('batch-tissue');
    const opt0 = document.getElementById('batch-tissue-opt0');
    tSel.disabled = !isHuman;
    if (!isHuman) tSel.value = '';
    if (opt0) opt0.textContent = isHuman ? 'Any / Pan-Tissue (Peak GTEx)' : 'N/A (Human GTEx only)';
    document.getElementById('batch-genes').placeholder = meta.batchPlaceholder;
    if (clearInputs) setCellValue('batch', '');
  } else {
    const tSel = document.getElementById('tissue-input');
    const opt0 = document.getElementById('single-tissue-opt0');
    tSel.disabled = !isHuman;
    if (!isHuman) tSel.value = '';
    if (opt0) opt0.textContent = isHuman ? 'Any / Pan-Tissue (Peak GTEx)' : 'N/A (Human GTEx only)';
    document.getElementById('gene-input').placeholder = meta.placeholder;
    if (clearInputs) {
      document.getElementById('gene-input').value = '';
      setCellValue('single', '');
    }
    renderExampleChips();
  }
  closeAllCellMenus();
}

function getCellEls(scope) {
  return scope === 'batch'
    ? {
        input: document.getElementById('batch-cell'),
        menu: document.getElementById('batch-cell-menu'),
        clearBtn: document.getElementById('batch-cell-clear')
      }
    : {
        input: document.getElementById('cell-input'),
        menu: document.getElementById('single-cell-menu'),
        clearBtn: document.getElementById('single-cell-clear')
      };
}

function setCellValue(scope, val) {
  const { input, clearBtn } = getCellEls(scope);
  input.value = val || '';
  clearBtn.style.display = input.value.trim() ? 'inline-flex' : 'none';
}

function clearCellInput(scope) {
  setCellValue(scope, '');
  openCellMenu(scope, true);
  getCellEls(scope).input.focus();
}

function closeAllCellMenus() {
  document.getElementById('single-cell-menu').classList.remove('open');
  document.getElementById('batch-cell-menu').classList.remove('open');
}

function toggleCellMenu(ev, scope) {
  if (ev) ev.stopPropagation();
  const { menu, input } = getCellEls(scope);
  if (menu.classList.contains('open')) {
    menu.classList.remove('open');
  } else {
    openCellMenu(scope, true);
    input.focus();
  }
}

function onCellInput(scope) {
  const { input, clearBtn } = getCellEls(scope);
  clearBtn.style.display = input.value.trim() ? 'inline-flex' : 'none';
  openCellMenu(scope, false);
}

function onCellKeyDown(ev, scope) {
  if (ev.key === 'Escape') {
    getCellEls(scope).menu.classList.remove('open');
  }
}

function pickCellByIndex(scope, idx) {
  const catalog = getCatalogForScope(scope);
  const chosen = idx < 0 ? '' : (catalog[idx] ? catalog[idx].name : '');
  setCellValue(scope, chosen);
  getCellEls(scope).menu.classList.remove('open');
}

function acceptTypedCell(scope) {
  getCellEls(scope).menu.classList.remove('open');
}

function openCellMenu(scope, showAll) {
  const other = scope === 'batch' ? 'single' : 'batch';
  getCellEls(other).menu.classList.remove('open');

  const sp = getScopeSpecies(scope);
  const meta = SPECIES_META[sp];
  const catalog = getCatalogForScope(scope);
  const { input, menu } = getCellEls(scope);
  const q = showAll ? '' : input.value.trim().toLowerCase();
  const curUpper = input.value.trim().toUpperCase();

  let html = `<button type="button" class="combo-opt ${!curUpper ? 'selected' : ''}" onmousedown="event.preventDefault(); pickCellByIndex('${scope}', -1)">
    <span><strong>Any / No Specific Model</strong></span>
    <span class="combo-opt-sub">${esc(meta.cellSub)}</span>
  </button>`;

  let lastGrp = '';
  let matchCount = 0;
  catalog.forEach((item, idx) => {
    if (q) {
      const hay = `${item.name} ${item.sub} ${item.group}`.toLowerCase();
      const compactHay = hay.replace(/[^a-z0-9]/g, '');
      const compactQ = q.replace(/[^a-z0-9]/g, '');
      if (!hay.includes(q) && (!compactQ || !compactHay.includes(compactQ))) {
        return;
      }
    }
    matchCount++;
    if (item.group !== lastGrp) {
      lastGrp = item.group;
      html += `<div class="combo-grp">${esc(item.group)}</div>`;
    }
    const isSel = item.name.toUpperCase() === curUpper;
    html += `<button type="button" class="combo-opt ${isSel ? 'selected' : ''}" onmousedown="event.preventDefault(); pickCellByIndex('${scope}', ${idx})">
      <span><strong>${esc(item.name)}</strong></span>
      <span class="combo-opt-sub">${esc(item.sub)}</span>
    </button>`;
  });

  if (q && matchCount === 0) {
    html += `<div class="combo-grp">${esc(meta.customHint)}</div>
      <button type="button" class="combo-opt selected" onmousedown="event.preventDefault(); acceptTypedCell('${scope}')">
        <span>Use <strong>"${esc(input.value.trim())}"</strong></span>
        <span class="combo-opt-sub">Press Enter to run</span>
      </button>`;
  }

  menu.innerHTML = html;
  menu.classList.add('open');
}

document.addEventListener('click', ev => {
  const sWrap = document.getElementById('single-combo-wrap');
  const bWrap = document.getElementById('batch-combo-wrap');
  if (sWrap && !sWrap.contains(ev.target)) {
    document.getElementById('single-cell-menu').classList.remove('open');
  }
  if (bWrap && !bWrap.contains(ev.target)) {
    document.getElementById('batch-cell-menu').classList.remove('open');
  }
});

function esc(s) {
  return String(s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

function csvCell(val) {
  const s = String(val ?? '');
  if (s.includes(',') || s.includes('"') || s.includes('\n')) {
    return '"' + s.replace(/"/g, '""') + '"';
  }
  return s;
}

function triggerCsvDownload(filename, lines) {
  const blob = new Blob([lines.join('\r\n') + '\r\n'], { type: 'text/csv;charset=utf-8' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
}

function switchAppMode(mode) {
  APP_MODE = mode;
  const isSingle = mode === 'SINGLE';
  document.getElementById('single-mode-view').style.display = isSingle ? 'block' : 'none';
  document.getElementById('batch-mode-view').style.display = isSingle ? 'none' : 'block';
  document.getElementById('tab-single').classList.toggle('active', isSingle);
  document.getElementById('tab-batch').classList.toggle('active', !isSingle);

  if (isSingle) {
    document.getElementById('batch-report-nav').style.display = 'none';
    if (SINGLE_PAYLOAD) {
      applyLoadedPayload(SINGLE_PAYLOAD);
    } else {
      document.getElementById('results-area').style.display = 'none';
    }
  } else {
    const selRow = BATCH_ROWS[SELECTED_BATCH_IDX];
    if (selRow && selRow.payload) {
      applyLoadedPayload(selRow.payload);
      renderBatchReportNav();
    } else {
      document.getElementById('results-area').style.display = 'none';
    }
  }
}

function flashCopied(btn, msg) {
  if (!btn) return;
  const old = btn.dataset.origLabel || btn.textContent;
  btn.dataset.origLabel = old;
  btn.textContent = msg || 'Copied!';
  setTimeout(() => { btn.textContent = old; }, 1300);
}

function copyText(txt, btn, msg) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(txt).then(() => flashCopied(btn, msg)).catch(() => fallbackCopyText(txt, btn, msg));
  } else {
    fallbackCopyText(txt, btn, msg);
  }
}

function fallbackCopyText(txt, btn, msg) {
  const ta = document.createElement('textarea');
  ta.value = txt;
  ta.style.position = 'fixed';
  ta.style.opacity = '0';
  document.body.appendChild(ta);
  ta.select();
  try { document.execCommand('copy'); } catch (e) {}
  document.body.removeChild(ta);
  flashCopied(btn, msg);
}

function syncUrlPermalink(gene, tissue, cell, species) {
  try {
    const params = new URLSearchParams();
    if (species && species !== 'human') params.set('species', species);
    if (gene) params.set('gene', gene);
    if (tissue) params.set('tissue', tissue);
    if (cell) params.set('cell_line', cell);
    const nextUrl = `${window.location.pathname}?${params.toString()}`;
    window.history.replaceState(null, '', nextUrl);
  } catch (e) {}
}

function copyShareableLink(btn) {
  copyText(window.location.href, btn, 'Link Copied!');
}

function copyExonRegion(exonNum, btn) {
  if (!CURRENT) return;
  const ex = CURRENT.report.exons.find(e => e.exon_number === exonNum);
  if (ex && ex.crispick_target_region) {
    copyText(ex.crispick_target_region, btn);
  }
}

function copyGuideSpacer(idx, btn) {
  const g = getDisplayedGuides()[idx];
  if (g) copyText(g.spacer_20nt, btn);
}

function copyOligoPair(idx, btn) {
  const g = getDisplayedGuides()[idx];
  if (g) copyText(["Fwd: 5'-" + g.cloning_fwd_oligo + "-3'", "Rev: 5'-" + g.cloning_rev_oligo + "-3'"].join('\r\n'), btn);
}

function copyPrimerPair(idx, btn) {
  const g = getDisplayedGuides()[idx];
  if (!g || !g.genotyping_primers) return;
  const p = g.genotyping_primers;
  const lines = [
    `${g.guide_id}_PCR_Fwd (${p.fwd_tm.toFixed(1)}C): 5'-${p.fwd_seq}-3'`,
    `${g.guide_id}_PCR_Rev (${p.rev_tm.toFixed(1)}C): 5'-${p.rev_seq}-3'`,
    `Amplicon: ${p.amplicon_bp} bp (${p.genomic_region}) | Sequence with ${p.recommended_seq_primer} primer`
  ];
  copyText(lines.join('\r\n'), btn);
}

function pinGuideToBasket(idx, btn) {
  if (!CURRENT) return;
  const g = getDisplayedGuides()[idx];
  if (!g) return;
  const exists = ORDER_BASKET.some(item => item.guide.guide_id === g.guide_id && item.guide.spacer_20nt === g.spacer_20nt);
  if (!exists) {
    ORDER_BASKET.push({ gene: CURRENT.report.gene_symbol, guide: g });
  }
  updateBasketUi();
  flashCopied(btn, 'Pinned!');
}

function clearOrderBasket() {
  ORDER_BASKET = [];
  updateBasketUi();
  if (CURRENT) renderGuidesTable();
}

function updateBasketUi() {
  const btn = document.getElementById('basket-btn');
  const clr = document.getElementById('basket-clear-btn');
  if (!ORDER_BASKET.length) {
    btn.style.display = 'none';
    clr.style.display = 'none';
    return;
  }
  const genes = [...new Set(ORDER_BASKET.map(x => x.gene))];
  btn.style.display = 'inline-block';
  clr.style.display = 'inline-block';
  btn.textContent = `Export Basket (${ORDER_BASKET.length} guides: ${genes.join('+')})`;
}

function buildIdtRowsFromItems(items) {
  const rows = ['Name,Sequence,Scale,Purification'];
  const seenPrimers = new Set();
  items.forEach(({ guide }) => {
    if (guide.cloning_fwd_oligo && guide.cloning_rev_oligo) {
      rows.push([`${guide.guide_id}_U6_Fwd`, guide.cloning_fwd_oligo, '25nm', 'STD'].map(csvCell).join(','));
      rows.push([`${guide.guide_id}_U6_Rev`, guide.cloning_rev_oligo, '25nm', 'STD'].map(csvCell).join(','));
    }
    const pr = guide.genotyping_primers;
    if (pr && pr.fwd_seq && pr.rev_seq) {
      const pKey = `${pr.fwd_seq}_${pr.rev_seq}`;
      if (!seenPrimers.has(pKey)) {
        seenPrimers.add(pKey);
        const baseName = guide.exon_number > 0 ? `${guide.guide_id.split('_E')[0]}_Exon${guide.exon_number}_TIDE` : `${guide.guide_id}_TIDE`;
        rows.push([`${baseName}_PCR_Fwd`, pr.fwd_seq, '25nm', 'STD'].map(csvCell).join(','));
        rows.push([`${baseName}_PCR_Rev`, pr.rev_seq, '25nm', 'STD'].map(csvCell).join(','));
      }
    }
  });
  return rows;
}

function exportIdtOrderCsv() {
  if (!CURRENT) return;
  const sym = CURRENT.report.gene_symbol;
  const currentItems = getDisplayedGuides().map(g => ({ gene: sym, guide: g }));
  const merged = [...currentItems];
  ORDER_BASKET.forEach(b => {
    if (!merged.some(m => m.guide.guide_id === b.guide.guide_id && m.guide.spacer_20nt === b.guide.spacer_20nt)) {
      merged.push(b);
    }
  });
  if (!merged.length) return;
  triggerCsvDownload(`${sym}_IDT_oligo_and_primer_order.csv`, buildIdtRowsFromItems(merged));
}

function exportBasketOnlyCsv() {
  if (!ORDER_BASKET.length) return;
  const genes = [...new Set(ORDER_BASKET.map(x => x.gene))].join('_');
  triggerCsvDownload(`${genes}_dual_IDT_order.csv`, buildIdtRowsFromItems(ORDER_BASKET));
}

async function inspectParalogByIndex(idx) {
  if (!CURRENT) return;
  const p = (CURRENT.report.paralogs || [])[idx];
  if (!p) return;
  const curSpecies = CURRENT.report.species_id || 'human';
  if (APP_MODE === 'BATCH') {
    const existingIdx = BATCH_ROWS.findIndex(r => r.gene.toUpperCase() === p.gene_symbol.toUpperCase());
    if (existingIdx >= 0 && BATCH_ROWS[existingIdx].payload) {
      openBatchRow(existingIdx, true);
      return;
    }
    const tissue = document.getElementById('batch-tissue').value.trim();
    const cell = document.getElementById('batch-cell').value.trim();
    const newIdx = BATCH_ROWS.length;
    BATCH_ROWS.push({ gene: p.gene_symbol, payload: null, error: null, pending: true });
    renderBatchTable();
    try {
      const params = new URLSearchParams({ _route: 'check', species: curSpecies, gene: p.gene_symbol });
      if (tissue) params.set('tissue', tissue);
      if (cell) params.set('cell_line', cell);
      const resp = await fetch(`/api/check?${params.toString()}`, { cache: 'no-store' });
      const data = await resp.json();
      if (!resp.ok || data.error) {
        BATCH_ROWS[newIdx].error = data.error || `HTTP ${resp.status}`;
      } else {
        BATCH_ROWS[newIdx].payload = data;
      }
    } catch (e) {
      BATCH_ROWS[newIdx].error = e.message;
    } finally {
      BATCH_ROWS[newIdx].pending = false;
      renderBatchTable();
      if (BATCH_ROWS[newIdx].payload) {
        openBatchRow(newIdx, true);
      }
    }
    return;
  }
  const tissue = document.getElementById('tissue-input').value.trim();
  const cell = document.getElementById('cell-input').value.trim();
  quickRun(p.gene_symbol, tissue, cell, curSpecies);
  window.scrollTo({ top: 0, behavior: 'smooth' });
}

function quickRun(gene, tissue, cell, species) {
  switchAppMode('SINGLE');
  if (species && SPECIES_META[species]) {
    document.getElementById('species-input').value = species;
    onSpeciesChange('single', false);
  }
  const sp = getScopeSpecies('single');
  document.getElementById('gene-input').value = gene;
  document.getElementById('tissue-input').value = sp === 'human' ? (tissue || '') : '';
  setCellValue('single', cell || '');
  closeAllCellMenus();
  syncUrlPermalink(gene, sp === 'human' ? (tissue || '') : '', cell || '', sp);
  runCheck(null, false);
}

async function runCheck(ev, updateUrl = true) {
  if (ev) ev.preventDefault();
  closeAllCellMenus();
  const sp = getScopeSpecies('single');
  const meta = SPECIES_META[sp];
  const gene = document.getElementById('gene-input').value.trim();
  const tissue = sp === 'human' ? document.getElementById('tissue-input').value.trim() : '';
  const cell = document.getElementById('cell-input').value.trim();
  if (!gene) return;

  if (updateUrl) {
    syncUrlPermalink(gene, tissue, cell, sp);
  }
  const reqId = ++ACTIVE_REQ_ID;
  const btn = document.getElementById('run-btn');
  const sbox = document.getElementById('status-box');
  btn.disabled = true;
  btn.textContent = 'Checking...';
  sbox.style.display = 'block';
  sbox.className = 'status-box status-loading';
  const tStart = performance.now();
  const renderLoadingMsg = () => {
    const secs = ((performance.now() - tStart) / 1000).toFixed(1);
    const dbDetail = sp === 'human'
      ? 'Querying GRCh38 isoforms, Pfam domains, DepMap, CCLE, gnomAD, GTEx &amp; HPA in parallel.'
      : `Querying ${esc(meta.assembly)} RefSeq/Ensembl exons, Pfam domains, Alliance DIOPT paralogs &amp; conserved human ortholog metrics in parallel.`;
    sbox.innerHTML = `<strong>Running ${esc(meta.label)} flight check on ${esc(gene)} (${secs}s)...</strong> ${dbDetail}`;
  };
  renderLoadingMsg();
  const timerHandle = setInterval(renderLoadingMsg, 150);

  try {
    const params = new URLSearchParams({ _route: 'check', species: sp, gene });
    if (tissue) params.set('tissue', tissue);
    if (cell) params.set('cell_line', cell);
    const resp = await fetch(`/api/check?${params.toString()}`, { cache: 'no-store' });
    const payload = await resp.json();
    if (reqId !== ACTIVE_REQ_ID) return;
    if (!resp.ok || payload.error) {
      throw new Error(payload.error || `HTTP ${resp.status}`);
    }
    SINGLE_PAYLOAD = payload;
    applyLoadedPayload(payload);
    sbox.className = 'status-box';
    sbox.style.display = 'none';
  } catch (err) {
    if (reqId !== ACTIVE_REQ_ID) return;
    sbox.style.display = 'block';
    sbox.className = 'status-box status-error';
    sbox.innerHTML = `<strong>Could not complete flight check for "${esc(gene)}":</strong> ${esc(err.message)}`;
  } finally {
    clearInterval(timerHandle);
    if (reqId === ACTIVE_REQ_ID) {
      btn.disabled = false;
      btn.textContent = 'Run Flight Check';
    }
  }
}

function applyLoadedPayload(payload) {
  CURRENT = payload;
  const coding = CURRENT.report.exons.filter(e => e.is_coding);
  const best = [...coding].sort((a, b) => b.priority_score - a.priority_score || a.exon_number - b.exon_number)[0];
  SELECTED_EXON = best ? best.exon_number : (coding[0] && coding[0].exon_number);
  // Auto-select CRISPRi scope when CRISPRi is the top recommended perturbation and promoter guides exist
  const modTxt = (CURRENT.report.modality.primary_modality || '').toUpperCase();
  if (modTxt.includes('CRISPRI') && (CURRENT.report.crispri_guides || []).length > 0) {
    GUIDE_SCOPE = 'CRISPRI';
  } else {
    GUIDE_SCOPE = 'PANEL';
  }
  document.getElementById('results-area').style.display = 'block';
  renderReport();
  triggerAsyncSnpCheck();
}

async function triggerAsyncSnpCheck() {
  if (!CURRENT) return;
  const snpReq = ++ACTIVE_SNP_REQ;
  const r = CURRENT.report;
  if ((r.species_id || 'human') !== 'human') return;
  const allGuides = [];
  (r.recommended_guides || []).forEach(g => allGuides.push(g));
  (r.crispri_guides || []).forEach(g => allGuides.push(g));
  (r.exons || []).forEach(ex => {
    (ex.candidate_guides || []).forEach(g => allGuides.push(g));
  });

  const specs = [];
  const seenIds = new Set();
  allGuides.forEach(g => {
    if (!g.guide_id || seenIds.has(g.guide_id) || SNP_CACHE[g.guide_id]) return;
    if (g.genomic_start_23bp && g.genomic_end_23bp) {
      seenIds.add(g.guide_id);
      specs.push({
        guide_id: g.guide_id,
        chrom: g.chrom,
        start: g.genomic_start_23bp,
        end: g.genomic_end_23bp,
        strand: g.genomic_strand
      });
    }
  });
  if (!specs.length) return;

  // Query in batches of up to 20 guides non-blockingly after report render
  try {
    const resp = await fetch('/api/snp-check?_route=snp-check', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ specs: specs.slice(0, 20) })
    });
    if (!resp.ok) return;
    const data = await resp.json();
    if (snpReq !== ACTIVE_SNP_REQ) return;
    if (data && data.results) {
      Object.assign(SNP_CACHE, data.results);
      renderGuidesTable();
    }
  } catch (e) {}
}

function selectExon(num) {
  SELECTED_EXON = num;
  renderSvgAndInspector();
  renderExonTable();
  if (GUIDE_SCOPE === 'EXON') {
    renderGuidesTable();
  } else {
    updateGuideScopeButtons();
  }
}

function selectDomain(aaStart, aaEnd) {
  if (!CURRENT) return;
  const coding = CURRENT.report.exons.filter(e => e.is_coding);
  const overlapping = coding.filter(e => Math.max(e.aa_start || 1, aaStart) <= Math.min(e.aa_end || 1, aaEnd));
  if (overlapping.length) {
    overlapping.sort((a, b) => b.priority_score - a.priority_score || a.exon_number - b.exon_number);
    selectExon(overlapping[0].exon_number);
  }
}

function toggleExonFilter() {
  ONLY_REC = !ONLY_REC;
  document.getElementById('filter-btn').textContent = ONLY_REC ? 'Showing: Recommended Only' : 'Showing: All Coding Exons';
  renderExonTable();
}

function setGuideScope(scope) {
  GUIDE_SCOPE = scope;
  renderGuidesTable();
}

function updateGuideScopeButtons() {
  if (!CURRENT) return;
  const r = CURRENT.report;
  const ciCount = (r.crispri_guides || []).length;
  const modRec = (r.modality.primary_modality || '').toUpperCase().includes('CRISPRI');
  document.getElementById('scope-panel').classList.toggle('active', GUIDE_SCOPE === 'PANEL');
  document.getElementById('scope-exon').classList.toggle('active', GUIDE_SCOPE === 'EXON');
  document.getElementById('scope-exon').textContent = SELECTED_EXON ? `Exon ${SELECTED_EXON} KO` : 'Selected Exon KO';
  const ciBtn = document.getElementById('scope-crispri');
  ciBtn.classList.toggle('active', GUIDE_SCOPE === 'CRISPRI');
  ciBtn.textContent = `CRISPRi Promoter (${ciCount})${modRec ? ' ★' : ''}`;
}

function getDisplayedGuides() {
  if (!CURRENT) return [];
  const r = CURRENT.report;
  if (GUIDE_SCOPE === 'CRISPRI') {
    return r.crispri_guides || [];
  }
  if (GUIDE_SCOPE === 'EXON') {
    const ex = r.exons.find(e => e.exon_number === SELECTED_EXON);
    return (ex && ex.candidate_guides) ? ex.candidate_guides : [];
  }
  return r.recommended_guides || [];
}

function copyGuidesFasta(btn) {
  const guides = getDisplayedGuides();
  if (!guides.length) return;
  const fasta = guides.map(g => {
    const loc = g.mode === 'CRISPRi'
      ? `TSS${g.tss_offset_bp >= 0 ? '+' : ''}${g.tss_offset_bp}bp`
      : `Exon${g.exon_number}_aa${g.aa_cut}`;
    return `>${g.guide_id} | ${loc} ${g.chrom}:${g.genomic_cut_pos} (${g.genomic_strand}) PAM=${g.pam}\r\n${g.spacer_20nt}`;
  }).join('\r\n');
  copyText(fasta, btn);
}

function renderSvgAndInspector() {
  const r = CURRENT.report;
  const coding = r.exons.filter(e => e.is_coding);
  const protLen = Math.max(1, r.primary_transcript.protein_length_aa);
  const W = 900, padX = 42, usable = W - 2 * padX;
  const aaToX = aa => padX + ((Math.max(1, Math.min(protLen, aa)) - 1) / Math.max(1, protLen - 1)) * usable;
  const fills = { RECOMMENDED: '#10b981', CAUTION: '#f59e0b', AVOID: '#f43f5e' };

  let exSvg = '';
  coding.forEach(ex => {
    const x1 = aaToX(ex.aa_start || 1);
    const x2 = Math.max(x1 + 7, aaToX(ex.aa_end || 1));
    const isSel = ex.exon_number === SELECTED_EXON;
    const stroke = isSel ? '#0f172a' : '#ffffff';
    const sw = isSel ? 2.5 : 1;
    const lbl = (x2 - x1) >= 24 ? `<text x="${((x1 + x2) / 2).toFixed(1)}" y="52" text-anchor="middle" font-size="12" font-weight="700" fill="#ffffff" pointer-events="none">E${ex.exon_number}</text>` : '';
    exSvg += `<g style="cursor:pointer" onclick="selectExon(${ex.exon_number})">
      <title>Exon ${ex.exon_number}: ${ex.priority_score}/100 (${ex.tier}) | aa ${ex.aa_start}-${ex.aa_end}</title>
      <rect x="${x1.toFixed(1)}" y="32" width="${(x2 - x1).toFixed(1)}" height="30" rx="4" fill="${fills[ex.tier] || '#94a3b8'}" stroke="${stroke}" stroke-width="${sw}"/>
      ${lbl}
    </g>`;
  });

  const doms = r.domains || [];
  const laneEnds = [];
  let domSvg = '';
  doms.forEach((dm, idx) => {
    const dNum = idx + 1;
    const x1 = aaToX(dm.aa_start);
    const x2 = Math.max(x1 + 18, aaToX(dm.aa_end));
    let lane = laneEnds.findIndex(endX => x1 >= endX + 4);
    if (lane === -1) {
      lane = laneEnds.length;
      laneEnds.push(x2);
    } else {
      laneEnds[lane] = x2;
    }
    const y = 95 + lane * 28;
    const w = x2 - x1;
    const fullText = `D${dNum}: ${dm.description}`;
    let barText = '';
    if (w >= fullText.length * 6.8 + 14) {
      barText = fullText;
    } else if (w >= 22) {
      barText = `D${dNum}`;
    }
    const lbl = barText ? `<text x="${((x1 + x2) / 2).toFixed(1)}" y="${y + 16}" text-anchor="middle" font-size="11.5" font-weight="700" fill="#ffffff" pointer-events="none">${esc(barText)}</text>` : '';
    domSvg += `<g style="cursor:pointer" onclick="selectDomain(${dm.aa_start}, ${dm.aa_end})">
      <title>D${dNum}: ${esc(dm.description)} [${esc(dm.source_db)}:${esc(dm.accession)}] (aa ${dm.aa_start}–${dm.aa_end}) — click to select overlapping exon</title>
      <rect x="${x1.toFixed(1)}" y="${y}" width="${w.toFixed(1)}" height="24" rx="4" fill="#4f46e5" fill-opacity="0.92" stroke="#312e81" stroke-width="1"/>
      ${lbl}
    </g>`;
  });

  const numDomLanes = Math.max(1, laneEnds.length);
  const domTrackHeight = numDomLanes * 28 - 4;
  if (!doms.length) {
    domSvg = `<text x="${(W / 2).toFixed(1)}" y="111" text-anchor="middle" font-size="12" fill="#64748b">No curated Pfam/Smart structural domains annotated on this isoform</text>`;
  }

  const track3LabelY = 95 + domTrackHeight + 23;
  const track3RectY = track3LabelY + 10;
  const axisY = track3RectY + 34;
  const svgH = axisY + 10;

  let epiSvg = '';
  (r.validation.epitope_windows || []).forEach(w => {
    const x1 = aaToX(w.aa_start);
    const x2 = Math.max(x1 + 5, aaToX(w.aa_end));
    epiSvg += `<g>
      <title>${w.status} (aa ${w.aa_start}-${w.aa_end}): ${esc(w.reason)}</title>
      <rect x="${x1.toFixed(1)}" y="${track3RectY}" width="${(x2 - x1).toFixed(1)}" height="16" rx="3" fill="${w.status === 'SAFE' ? '#0d9488' : '#f43f5e'}" fill-opacity="0.85"/>
    </g>`;
  });

  document.getElementById('svg-wrap').innerHTML = `
    <svg viewBox="0 0 ${W} ${svgH}" style="width:100%;min-width:640px;">
      <text x="${padX}" y="22" font-size="12" font-weight="700" fill="#334155">Track 1: Coding Exons</text>
      <line x1="${padX}" y1="47" x2="${W - padX}" y2="47" stroke="#cbd5e1" stroke-width="3"/>
      ${exSvg}
      <text x="${padX}" y="85" font-size="12" font-weight="700" fill="#334155">Track 2: Conserved Pfam / Structural Domains (deduplicated; click any domain to jump to its exon)</text>
      <rect x="${padX}" y="95" width="${usable}" height="${domTrackHeight}" rx="4" fill="#f1f5f9"/>
      ${domSvg}
      <text x="${padX}" y="${track3LabelY}" font-size="12" font-weight="700" fill="#334155">Track 3: Western Antibody Epitope Map (Teal = Safe | Coral = Re-initiation / Skippable Risk)</text>
      <rect x="${padX}" y="${track3RectY}" width="${usable}" height="16" rx="3" fill="#f1f5f9"/>
      ${epiSvg}
      <text x="${padX}" y="${axisY}" font-size="12" fill="#64748b">N-term (aa 1)</text>
      <text x="${(W / 2).toFixed(1)}" y="${axisY}" text-anchor="middle" font-size="12" fill="#64748b">aa ${Math.floor(protLen / 2)}</text>
      <text x="${W - padX}" y="${axisY}" text-anchor="end" font-size="12" fill="#64748b">C-term (aa ${protLen})</text>
    </svg>`;

  const domLegendEl = document.getElementById('domain-legend');
  if (doms.length) {
    domLegendEl.innerHTML = `<div class="dom-legend">
      <span style="font-size:12px;font-weight:700;color:#312e81;">Track 2 Domain Key (full names):</span>
      ${doms.map((dm, idx) => `
        <button type="button" class="dom-chip" onclick="selectDomain(${dm.aa_start}, ${dm.aa_end})" title="Click to inspect the best coding exon overlapping ${esc(dm.description)}">
          <span class="dom-badge">D${idx + 1}</span>
          <strong>${esc(dm.description)}</strong>
          <span style="color:#475569;">(${esc(dm.source_db)}:${esc(dm.accession)} &bull; aa ${dm.aa_start}–${dm.aa_end})</span>
        </button>
      `).join('')}
    </div>`;
  } else {
    domLegendEl.innerHTML = '';
  }

  const ex = coding.find(e => e.exon_number === SELECTED_EXON) || coding[0];
  if (ex) {
    document.getElementById('exon-inspector').innerHTML = `
      <div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:10px;">
        <div>
          <strong>Selected: Exon ${ex.exon_number}</strong>
          <span class="pill pill-${ex.tier}" style="margin-left:8px;">${ex.priority_score}/100 &bull; ${ex.tier}</span>
          <span style="font-size:12px;color:var(--muted);margin-left:8px;">aa ${ex.aa_start}–${ex.aa_end} &bull; ${ex.exon_length_bp} bp (${ex.is_symmetric ? 'divisible by 3 — symmetric' : 'asymmetric'}) &bull; ${ex.spcas9_pam_count ?? '—'} SpCas9 NGG cut sites</span>
        </div>
        <div style="display:flex;align-items:center;gap:8px;">
          <code>${esc(ex.crispick_target_region)}</code>
          <button type="button" class="btn-secondary" onclick="copyExonRegion(${ex.exon_number}, this)">Copy Region</button>
          <button type="button" class="btn-secondary" onclick="setGuideScope('EXON')">View Exon ${ex.exon_number} Guides</button>
        </div>
      </div>
      <ul style="margin:8px 0 0 0;padding-left:20px;font-size:12.5px;color:#334155;">
        ${ex.score_breakdown.map(b => `<li>${esc(b)}</li>`).join('')}
      </ul>`;
  }
}

function renderSnpBadge(guideId) {
  if (CURRENT && (CURRENT.report.species_id || 'human') !== 'human') {
    return '<span class="pill pill-INFO" title="Sequence your wild-type lab strain or stock across the target exon with the paired PCR primers">Seq WT Strain</span>';
  }
  const info = SNP_CACHE[guideId];
  if (!info) {
    return '<span class="pill pill-INFO" style="font-weight:600;">Checking gnomAD...</span>';
  }
  if (info.status === 'CLEAN') {
    return '<span class="pill pill-RECOMMENDED" title="Zero common human variants (MAF >= 1%) in 20-nt spacer + PAM">&#10003; Clean (&lt;1% MAF)</span>';
  }
  if (info.status === 'HIGH') {
    return `<span class="pill pill-AVOID" title="${esc(info.summary)}">&#9888; ${esc(info.summary)}</span>`;
  }
  if (info.status === 'WARN') {
    return `<span class="pill pill-CAUTION" title="${esc(info.summary)}">&#9888; ${esc(info.summary)}</span>`;
  }
  return '<span style="font-size:11.5px;color:var(--muted);">—</span>';
}

function renderGuidesTable() {
  updateGuideScopeButtons();
  const guides = getDisplayedGuides();
  const wrap = document.getElementById('guides-wrap');
  const isCrispri = GUIDE_SCOPE === 'CRISPRI';

  if (!guides.length) {
    const emptyMsg = isCrispri
      ? 'No canonical promoter SpCas9 guides designed in TSS [-50, +300 bp] for this locus.'
      : 'No clean non-<code>TTTT</code> canonical SpCas9 <code>NGG</code> spacers found in this scope (consider an SpG/SpRY PAM-relaxed Cas9 variant or select another exon).';
    wrap.innerHTML = `<div style="font-size:13px;color:var(--muted);padding:8px 0;">${emptyMsg}</div>`;
    return;
  }

  const tssBanner = (isCrispri && CURRENT.report.tss_genomic_pos)
    ? `<div style="background:#f0fdfa;border:1px solid #99f6e4;color:#0f766e;border-radius:8px;padding:8px 12px;font-size:12.5px;margin-bottom:10px;">
         <strong>dCas9-KRAB Promoter Window:</strong> Guides designed around primary RefSeq 5&prime; TSS at <code>${esc(CURRENT.report.chrom)}:${CURRENT.report.tss_genomic_pos}</code> (Horlbeck/Weissman optimal window: <code>-50</code> to <code>+300 bp</code>; peak <code>+25..+150 bp</code>). No DSBs — avoids multi-copy amplification &amp; pan-essential clonal escaper artifacts.
       </div>`
    : '';

  wrap.innerHTML = `${tssBanner}
    <table>
      <thead>
        <tr>
          <th>Guide ID</th>
          <th>${isCrispri ? 'TSS Offset' : 'Exon &amp; AA'}</th>
          <th>20-nt Spacer (5&prime;&rarr;3&prime;) + PAM</th>
          <th>Common SNP Shield</th>
          <th>GC% &amp; Quality</th>
          <th>U6 Golden Gate Oligos</th>
          <th>${isCrispri ? 'Order Basket' : 'TIDE / ICE PCR Primers &amp; Basket'}</th>
        </tr>
      </thead>
      <tbody>
        ${guides.map((g, idx) => {
          const domPill = g.overlaps_domain ? `<div style="font-size:11px;color:#4338ca;font-weight:600;">&#9670; ${esc(g.overlaps_domain)}</div>` : '';
          const locCell = isCrispri
            ? `<code>${g.tss_offset_bp >= 0 ? '+' : ''}${g.tss_offset_bp} bp</code><div style="font-size:11px;color:var(--muted);">${esc(g.chrom)}:${g.genomic_cut_pos} (${esc(g.genomic_strand)})</div>`
            : `Exon ${g.exon_number} &bull; aa ${g.aa_cut}${domPill}<div style="font-size:11px;color:var(--muted);">${esc(g.chrom)}:${g.genomic_cut_pos} (${esc(g.genomic_strand)})</div>`;
          const notesTxt = (g.quality_notes || []).slice(0, 2).join(' • ');
          const isPinned = ORDER_BASKET.some(b => b.guide.guide_id === g.guide_id && b.guide.spacer_20nt === g.spacer_20nt);
          const pinBtn = `<button type="button" class="btn-secondary" onclick="pinGuideToBasket(${idx}, this)">${isPinned ? '&#10003; Pinned' : '+ Pin Guide'}</button>`;

          let primerCell = pinBtn;
          if (!isCrispri && g.genotyping_primers) {
            const p = g.genotyping_primers;
            primerCell = `
              <div style="font-size:11.5px;white-space:nowrap;">
                <div><strong>${p.amplicon_bp} bp</strong> &bull; Seq w/ <strong>${esc(p.recommended_seq_primer)}</strong></div>
                <div style="color:var(--muted);">Fwd (${p.fwd_tm.toFixed(0)}&deg;C): <code>${esc(p.fwd_seq)}</code></div>
                <div style="color:var(--muted);">Rev (${p.rev_tm.toFixed(0)}&deg;C): <code>${esc(p.rev_seq)}</code></div>
                <div style="margin-top:4px;display:flex;gap:5px;">
                  <button type="button" class="btn-secondary" onclick="copyPrimerPair(${idx}, this)">Copy PCR</button>
                  ${pinBtn}
                </div>
              </div>`;
          }

          return `<tr>
            <td><strong>${esc(g.guide_id)}</strong></td>
            <td style="white-space:nowrap;">${locCell}</td>
            <td style="white-space:nowrap;">
              <code style="font-size:13px;font-weight:700;color:#0f172a;">${esc(g.spacer_20nt)}</code>
              <span class="pam-badge">${esc(g.pam)}</span>
              <div style="margin-top:4px;"><button type="button" class="btn-secondary" onclick="copyGuideSpacer(${idx}, this)">Copy 20-nt</button></div>
            </td>
            <td>${renderSnpBadge(g.guide_id)}</td>
            <td><strong>${g.gc_percent}% GC</strong> &bull; ${g.sequence_score}/100<div style="font-size:11.5px;color:var(--muted);">${esc(notesTxt)}</div></td>
            <td style="font-size:11.5px;white-space:nowrap;">
              <div>Fwd: <code>${esc(g.cloning_fwd_oligo)}</code></div>
              <div style="margin-top:2px;">Rev: <code>${esc(g.cloning_rev_oligo)}</code></div>
              <button type="button" class="btn-secondary" style="margin-top:4px;" onclick="copyOligoPair(${idx}, this)">Copy Fwd+Rev</button>
            </td>
            <td>${primerCell}</td>
          </tr>`;
        }).join('')}
      </tbody>
    </table>`;
}

function renderExonTable() {
  const coding = CURRENT.report.exons.filter(e => e.is_coding);
  const rows = ONLY_REC ? coding.filter(e => e.tier === 'RECOMMENDED') : coding;
  document.getElementById('exon-tbody').innerHTML = rows.map(ex => `
    <tr class="clickable-row ${ex.exon_number === SELECTED_EXON ? 'selected-row' : ''}" onclick="selectExon(${ex.exon_number})">
      <td><strong>Exon ${ex.exon_number}</strong></td>
      <td><strong>${ex.priority_score}</strong>/100</td>
      <td><span class="pill pill-${ex.tier}">${ex.tier}</span></td>
      <td>aa ${ex.aa_start}–${ex.aa_end}</td>
      <td>${ex.exon_length_bp} bp ${ex.is_symmetric ? '(div3)' : '(asym)'}</td>
      <td>${Math.round(ex.isoform_coverage_fraction * 100)}%</td>
      <td>${ex.spcas9_pam_count ?? '—'}</td>
      <td><code>${esc(ex.nmd_status)}</code></td>
      <td><code>${esc(ex.crispick_target_region)}</code> <button type="button" class="btn-secondary" onclick="event.stopPropagation();copyExonRegion(${ex.exon_number}, this)">Copy</button></td>
    </tr>`).join('');
}

function renderReport() {
  const r = CURRENT.report;
  const tx = r.primary_transcript;
  const isHuman = (r.species_id || 'human') === 'human';
  const spCommon = r.species_common_name || 'Human';
  const spSci = r.species_scientific_name || 'Homo sapiens';
  const orthNote = (!isHuman && r.human_ortholog_symbol) ? ` • 1:1 Human Ortholog: ${r.human_ortholog_symbol}` : '';

  document.getElementById('res-title').innerHTML = `${esc(r.gene_symbol)} <span style="font-weight:400;color:var(--muted);font-size:16px;">(${esc(r.ensembl_id)})</span>`;
  document.getElementById('res-sub').textContent = `${spSci} (${spCommon})${orthNote} • ${r.gene_description} • ${r.assembly} ${r.chrom} (${r.strand >= 0 ? '+' : '-'} strand) • Primary transcript: ${tx.transcript_id} (${tx.protein_length_aa} aa)`;
  updateBasketUi();

  document.getElementById('lbl-kpi-1').textContent = isHuman
    ? 'DepMap Fitness & CCLE Genetics'
    : `Conserved Fitness${r.human_ortholog_symbol ? ` (via ${r.human_ortholog_symbol})` : ''}`;
  document.getElementById('lbl-kpi-2').textContent = isHuman
    ? 'Human LoF Constraint (gnomAD)'
    : 'Ortholog LoF Constraint (gnomAD)';
  document.getElementById('lbl-kpi-3').textContent = isHuman
    ? 'GTEx Baseline Expression'
    : 'Ortholog Baseline Expression';
  document.getElementById('exon-coord-th').textContent = `Target Region (${r.assembly})`;
  document.getElementById('paralog-hdr').textContent = `5. ${spCommon} Paralog Buffering Check`;

  const guidesSub = document.getElementById('guides-sub');
  if (isHuman) {
    guidesSub.innerHTML = 'Includes U6 Golden Gate oligos (<code>BsmBI</code>/<code>BbsI</code>), paired 430–740 bp TIDE/ICE PCR primers, and background gnomAD common SNP screening (MAF &ge;1%).';
  } else if (r.species_id === 'fly') {
    guidesSub.innerHTML = 'Includes U6 Golden Gate oligos (note: swap <code>5&prime;-CACC</code> &rarr; <code>5&prime;-GTCG</code> if cloning into Drosophila <code>pCFD3</code> <code>dU6:3</code>) and paired 430–740 bp TIDE/ICE PCR primers.';
  } else {
    guidesSub.innerHTML = `Includes U6 Golden Gate oligos (<code>BsmBI</code>/<code>BbsI</code>) and paired 430–740 bp TIDE/ICE PCR primers on <code>${esc(r.assembly)}</code>.`;
  }

  document.getElementById('offtarget-assembly-note').innerHTML = `
    <strong>Before ordering oligos — 10-second genome-wide off-target check:</strong> Local rules rank on-target quality &amp; biological safety; paste your chosen 2–3 spacers into <strong>Cas-OFFinder</strong>${(r.species_id === 'human' || r.species_id === 'mouse' || r.species_id === 'rat') ? ' or <strong>Broad CRISPick</strong>' : ''} to verify zero 0–2 mismatch off-target sites across <strong>${esc(r.assembly)}</strong>.`;

  const vEl = document.getElementById('res-verdict');
  vEl.className = `verdict verdict-${r.modality.risk_level}`;
  vEl.innerHTML = `
    <div style="font-size:12px;font-weight:800;text-transform:uppercase;">Recommended Perturbation Strategy (${esc(r.modality.risk_level)} Trap Risk)</div>
    <div style="font-size:20px;font-weight:800;margin:4px 0;">${esc(r.modality.primary_modality)}</div>
    <div style="font-size:14px;">${esc(r.modality.summary_verdict)}</div>`;

  const ess = r.essentiality;
  if (ess.matched_cell_line && ess.matched_cell_line_chronos !== null) {
    const c = ess.matched_cell_line_chronos;
    document.getElementById('k-chronos').textContent = `${ess.matched_cell_line}: ${c >= 0 ? '+' : ''}${c.toFixed(2)}`;
  } else if (ess.depmap_median_chronos !== null) {
    const c = ess.depmap_median_chronos;
    document.getElementById('k-chronos').textContent = `Median: ${c >= 0 ? '+' : ''}${c.toFixed(2)}`;
  } else {
    document.getElementById('k-chronos').textContent = '—';
  }
  document.getElementById('k-chronos-sub').textContent = ess.depmap_pan_essential_fraction !== null
    ? `${ess.depmap_essential_lines}/${ess.depmap_total_lines} lines (${Math.round(ess.depmap_pan_essential_fraction * 100)}%) <= -1.0`
    : 'No DepMap CRISPR screen data';

  // Render CCLE endogenous mutation / copy-number status badge
  const ccleEl = document.getElementById('k-ccle-sub');
  const clg = r.cell_line_genetics;
  if (clg && clg.matched_ccle_sample) {
    const shortCell = clg.matched_ccle_sample.split('_')[0];
    if (clg.has_truncating_lof || clg.cna_gistic === -2) {
      const detail = clg.mutations[0] || clg.cna_label || 'LoF';
      ccleEl.innerHTML = `<span class="ccle-tag ccle-alert">CCLE ${esc(shortCell)}: ${esc(detail)}</span>`;
    } else if (clg.cna_gistic === 2) {
      ccleEl.innerHTML = `<span class="ccle-tag ccle-alert">CCLE ${esc(shortCell)}: High Amplification (GISTIC +2)</span>`;
    } else if (clg.mutations && clg.mutations.length) {
      ccleEl.innerHTML = `<span class="ccle-tag ccle-warn">CCLE ${esc(shortCell)}: ${esc(clg.mutations[0])}</span>`;
    } else {
      const cnTxt = clg.cna_gistic === -1 ? 'Hemizygous (-1)' : (clg.cna_gistic === 1 ? 'Copy Gain (+1)' : 'Copy-Neutral');
      ccleEl.innerHTML = `<span class="ccle-tag ccle-ok">CCLE ${esc(shortCell)}: Coding WT &bull; ${cnTxt}</span>`;
    }
  } else {
    ccleEl.innerHTML = '';
  }

  document.getElementById('k-loeuf').textContent = ess.gnomad_loeuf !== null ? `LOEUF ${ess.gnomad_loeuf.toFixed(2)}` : '—';
  document.getElementById('k-pli-sub').textContent = ess.gnomad_pli !== null
    ? `gnomAD pLI = ${ess.gnomad_pli.toFixed(2)}`
    : (ess.gnomad_oe_lof !== null ? `o/e LoF = ${ess.gnomad_oe_lof.toFixed(2)} (<0.35 = constrained)` : '< 0.35 indicates strong constraint');

  const exp = r.expression;
  if (exp.resolved_gtex_tissue && exp.requested_tissue_tpm !== null) {
    document.getElementById('k-tpm').textContent = `${exp.requested_tissue_tpm} TPM`;
    document.getElementById('k-loc').textContent = `${exp.resolved_gtex_tissue} • ${(exp.hpa_subcellular_locations || []).slice(0, 2).join(', ') || 'HPA unannotated'}`;
  } else if (exp.top_tissues && exp.top_tissues.length) {
    document.getElementById('k-tpm').textContent = `${exp.top_tissues[0][1]} TPM (peak)`;
    document.getElementById('k-loc').textContent = `Top: ${exp.top_tissues[0][0]}`;
  } else {
    document.getElementById('k-tpm').textContent = '—';
    document.getElementById('k-loc').textContent = '';
  }

  document.getElementById('k-iso').textContent = `${1 + r.alternative_isoforms.length} isoform(s)`;
  document.getElementById('k-exons').textContent = `${tx.num_exons} total exon(s) on primary transcript`;

  renderSvgAndInspector();

  const warns = r.modality.warnings || [];
  document.getElementById('traps-wrap').innerHTML = warns.length ? warns.map(w => `
    <div class="trap trap-${w.severity}">
      <div style="font-weight:700;margin-bottom:4px;">[${esc(w.severity)}] ${esc(w.title)}</div>
      <div><strong>Why:</strong> ${esc(w.explanation)}</div>
      <div style="margin-top:3px;"><strong>Bench Action:</strong> ${esc(w.bench_action)}</div>
    </div>`).join('') : '<div style="font-size:13px;color:var(--muted);">No major biological KO failure traps flagged.</div>';

  renderExonTable();
  renderGuidesTable();

  const pars = r.paralogs || [];
  document.getElementById('paralogs-wrap').innerHTML = pars.length ? pars.map((p, idx) => `
    <div style="display:flex;justify-content:space-between;align-items:center;gap:8px;padding:9px 0;border-bottom:1px solid var(--border);font-size:13px;">
      <div>
        <strong>${esc(p.gene_symbol)}</strong> <code>${esc(p.ensembl_id)}</code>
        <div style="font-size:12px;color:var(--muted);">${p.percent_identity}% query / ${p.target_percent_identity}% target id${p.target_tissue_tpm !== null ? ` • Tissue: ${p.target_tissue_tpm} TPM` : ''}</div>
      </div>
      <div style="display:flex;align-items:center;gap:6px;white-space:nowrap;">
        <span class="pill pill-${p.buffering_risk === 'HIGH' ? 'AVOID' : (p.buffering_risk === 'MODERATE' ? 'CAUTION' : 'RECOMMENDED')}">${esc(p.buffering_risk)}</span>
        <button type="button" class="btn-secondary" onclick="inspectParalogByIndex(${idx})">Inspect &rarr;</button>
      </div>
    </div>`).join('') : `<p style="font-size:13px;color:var(--muted);">No ${esc(spCommon.toLowerCase())} paralogs >= 35% sequence identity.</p>`;

  const val = r.validation;
  const tools = val.orthogonal_chemical_tools || [];
  const toolsList = tools.length ? `<ul style="padding-left:18px;margin:6px 0 0 0;font-size:13px;">${tools.map(t => `
    <li style="margin-bottom:4px;"><strong>${esc(t.name)}</strong> (${esc(t.tool_type)}, ${esc(t.quality_or_stage)}) — ${esc(t.mechanism_of_action)} ${t.has_negative_control ? '<span class="pill pill-RECOMMENDED">+Neg Control</span>' : ''}</li>`).join('')}</ul>`
    : '<div style="font-size:13px;color:var(--muted);margin-top:4px;">None curated in Open Targets — rely on cDNA rescue + &ge;2 independent sgRNAs.</div>';

  document.getElementById('val-wrap').innerHTML = `
    <div style="font-size:13px;margin-bottom:8px;"><strong>Protein Readout:</strong> ${esc(val.recommended_protein_assay)}</div>
    <div style="font-size:13px;margin-bottom:8px;"><strong>Western Antibody:</strong> ${esc(val.epitope_guidance_summary)}</div>
    <div style="font-size:13px;margin-bottom:10px;"><strong>RT-PCR Splicing Check:</strong> ${esc(val.rt_qpcr_guidance)}</div>
    <div style="font-size:13px;font-weight:700;">Orthogonal Chemical Probes / Clinical Modulators:</div>
    ${toolsList}`;
}

function downloadReport(kind) {
  if (!CURRENT) return;
  const sym = CURRENT.report.gene_symbol;
  let content = '', mime = 'text/plain', ext = 'txt';
  if (kind === 'html') { content = CURRENT.standalone_html; mime = 'text/html'; ext = 'html'; }
  else if (kind === 'md') { content = CURRENT.markdown; mime = 'text/markdown'; ext = 'md'; }
  else { content = JSON.stringify(CURRENT.report, null, 2); mime = 'application/json'; ext = 'json'; }

  const blob = new Blob([content], { type: mime });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = `${sym}_crispr_preflight.${ext}`;
  a.click();
}

// ==================== BATCH HIT-LIST TRIAGE ====================

function loadDemoBatch() {
  const sp = getScopeSpecies('batch');
  const meta = SPECIES_META[sp];
  document.getElementById('batch-genes').value = meta.demoBatch;
  document.getElementById('batch-tissue').value = '';
  setCellValue('batch', '');
  runBatchTriage();
}

async function runBatchTriage() {
  const sp = getScopeSpecies('batch');
  const raw = document.getElementById('batch-genes').value || '';
  const seenLower = new Set();
  const tokens = [];
  raw.split(/[\s,;]+/).forEach(tok => {
    const clean = tok.trim();
    if (!clean) return;
    const key = clean.toLowerCase();
    if (seenLower.has(key)) return;
    seenLower.add(key);
    tokens.push(sp === 'human' ? clean.toUpperCase() : clean);
  });
  tokens.splice(30);
  if (!tokens.length) return;

  const tissue = sp === 'human' ? document.getElementById('batch-tissue').value.trim() : '';
  const cell = document.getElementById('batch-cell').value.trim();
  const btn = document.getElementById('batch-run-btn');
  const prog = document.getElementById('batch-progress');
  btn.disabled = true;
  btn.textContent = 'Triaging...';
  prog.style.display = 'block';
  document.getElementById('batch-results-card').style.display = 'block';

  BATCH_ROWS = tokens.map(g => ({ gene: g, payload: null, error: null, pending: true }));
  SELECTED_BATCH_IDX = -1;
  document.getElementById('results-area').style.display = 'none';
  renderBatchTable();

  let done = 0;
  let cursor = 0;
  const updateProg = () => {
    prog.textContent = done < tokens.length
      ? `Completed ${done}/${tokens.length} ${SPECIES_META[sp].label} genes (running 3 parallel checks)...`
      : `Finished all ${tokens.length} ${SPECIES_META[sp].label} genes! Click any gene row in the matrix below to view its full report right underneath.`;
  };
  updateProg();

  async function worker() {
    while (cursor < BATCH_ROWS.length) {
      const idx = cursor++;
      const row = BATCH_ROWS[idx];
      try {
        const params = new URLSearchParams({ _route: 'check', species: sp, gene: row.gene });
        if (tissue) params.set('tissue', tissue);
        if (cell) params.set('cell_line', cell);
        const resp = await fetch(`/api/check?${params.toString()}`, { cache: 'no-store' });
        const data = await resp.json();
        if (!resp.ok || data.error) {
          row.error = data.error || `HTTP ${resp.status}`;
        } else {
          row.payload = data;
        }
      } catch (e) {
        row.error = e.message;
      } finally {
        row.pending = false;
        done++;
        updateProg();
        if (APP_MODE === 'BATCH' && SELECTED_BATCH_IDX === -1 && row.payload) {
          SELECTED_BATCH_IDX = idx;
          applyLoadedPayload(row.payload);
        }
        renderBatchTable();
        if (APP_MODE === 'BATCH' && SELECTED_BATCH_IDX >= 0) {
          renderBatchReportNav();
        }
      }
    }
  }

  const concurrency = Math.min(3, BATCH_ROWS.length);
  await Promise.all(Array.from({ length: concurrency }, () => worker()));
  btn.disabled = false;
  btn.textContent = 'Triage Gene List';
}

function renderBatchTable() {
  const tbody = document.getElementById('batch-tbody');
  tbody.innerHTML = BATCH_ROWS.map((row, idx) => {
    if (row.pending) {
      return `<tr><td><strong>${esc(row.gene)}</strong></td><td colspan="8" style="color:var(--muted);">Running parallel flight check...</td></tr>`;
    }
    if (row.error || !row.payload) {
      return `<tr><td><strong>${esc(row.gene)}</strong></td><td colspan="8" style="color:#b91c1c;">Error: ${esc(row.error || 'Failed')}</td></tr>`;
    }
    const r = row.payload.report;
    const riskCls = r.modality.risk_level === 'HIGH' ? 'AVOID' : (r.modality.risk_level === 'MODERATE' ? 'CAUTION' : 'RECOMMENDED');
    const ess = r.essentiality;
    const chronosTxt = (ess.matched_cell_line && ess.matched_cell_line_chronos !== null)
      ? `${ess.matched_cell_line}: ${ess.matched_cell_line_chronos >= 0 ? '+' : ''}${ess.matched_cell_line_chronos.toFixed(2)}`
      : (ess.depmap_median_chronos !== null ? `Med: ${ess.depmap_median_chronos.toFixed(2)}` : '—');

    let ccleTxt = '—';
    const clg = r.cell_line_genetics;
    if (clg && clg.matched_ccle_sample) {
      if (clg.has_truncating_lof || clg.cna_gistic === -2) {
        ccleTxt = `<span style="color:#b91c1c;font-weight:700;">LoF: ${esc(clg.mutations[0] || clg.cna_label)}</span>`;
      } else if (clg.cna_gistic === 2) {
        ccleTxt = '<span style="color:#b91c1c;font-weight:700;">High Amp (+2)</span>';
      } else if (clg.mutations && clg.mutations.length) {
        ccleTxt = `<span style="color:#b45309;font-weight:600;">${esc(clg.mutations[0])}</span>`;
      } else {
        ccleTxt = 'WT / Neutral';
      }
    }

    const tpmTxt = r.expression.requested_tissue_tpm !== null
      ? `${r.expression.requested_tissue_tpm} TPM`
      : (r.expression.top_tissues && r.expression.top_tissues.length ? `${r.expression.top_tissues[0][1]} TPM (peak)` : '—');

    const topPar = (r.paralogs || [])[0];
    const parTxt = topPar ? `${esc(topPar.gene_symbol)} (${esc(topPar.buffering_risk)})` : 'None';

    const coding = r.exons.filter(e => e.is_coding).sort((a, b) => b.priority_score - a.priority_score);
    const bestEx = coding[0];
    const topG = (r.recommended_guides || [])[0] || (r.crispri_guides || [])[0];
    const guideTxt = bestEx
      ? `E${bestEx.exon_number} (${bestEx.priority_score}/100)${topG ? `<br/><code>${esc(topG.spacer_20nt)}</code>` : ''}`
      : '—';

    const isSel = idx === SELECTED_BATCH_IDX;
    const btnHtml = isSel
      ? `<button type="button" class="btn-primary" style="height:auto;padding:5px 10px;font-size:12px;" onclick="event.stopPropagation();openBatchRow(${idx}, true)">&#10003; Viewing Below</button>`
      : `<button type="button" class="btn-accent" onclick="event.stopPropagation();openBatchRow(${idx}, true)">View Below &darr;</button>`;

    return `<tr class="clickable-row ${isSel ? 'selected-row' : ''}" onclick="openBatchRow(${idx}, true)">
      <td><strong>${esc(r.gene_symbol)}</strong><div style="font-size:11px;color:var(--muted);">${esc(r.ensembl_id)}</div></td>
      <td><span class="pill pill-${riskCls}">${esc(r.modality.risk_level)}</span></td>
      <td><strong>${esc(r.modality.primary_modality)}</strong></td>
      <td style="white-space:nowrap;">${esc(chronosTxt)}</td>
      <td>${ccleTxt}</td>
      <td style="white-space:nowrap;">${esc(tpmTxt)}</td>
      <td>${parTxt}</td>
      <td>${guideTxt}</td>
      <td style="white-space:nowrap;">${btnHtml}</td>
    </tr>`;
  }).join('');
}

function renderBatchReportNav() {
  const nav = document.getElementById('batch-report-nav');
  if (APP_MODE !== 'BATCH' || SELECTED_BATCH_IDX < 0 || !BATCH_ROWS[SELECTED_BATCH_IDX] || !BATCH_ROWS[SELECTED_BATCH_IDX].payload) {
    nav.style.display = 'none';
    return;
  }
  const readyIndices = BATCH_ROWS.map((r, i) => (r.payload ? i : -1)).filter(i => i >= 0);
  const posInReady = readyIndices.indexOf(SELECTED_BATCH_IDX);
  const prevIdx = posInReady > 0 ? readyIndices[posInReady - 1] : -1;
  const nextIdx = (posInReady >= 0 && posInReady < readyIndices.length - 1) ? readyIndices[posInReady + 1] : -1;

  const chipsHtml = readyIndices.map(i => {
    const sym = BATCH_ROWS[i].payload.report.gene_symbol;
    const active = i === SELECTED_BATCH_IDX;
    const style = active
      ? 'background:#0f766e;color:#ffffff;border-color:#0f766e;'
      : 'background:#ffffff;color:#0f766e;border-color:#5eead4;';
    return `<button type="button" class="chip" style="${style}" onclick="openBatchRow(${i}, false)">${esc(sym)}</button>`;
  }).join('');

  nav.style.display = 'block';
  nav.innerHTML = `
    <div style="display:flex;flex-wrap:wrap;justify-content:space-between;align-items:center;gap:10px;">
      <div style="display:flex;flex-wrap:wrap;align-items:center;gap:6px;">
        <span style="font-size:12.5px;font-weight:800;color:#0f766e;margin-right:4px;">Detailed Report (${posInReady + 1}/${readyIndices.length}):</span>
        ${chipsHtml}
      </div>
      <div style="display:flex;gap:6px;align-items:center;">
        <button type="button" class="btn-secondary" ${prevIdx < 0 ? 'disabled style="opacity:0.45;"' : `onclick="openBatchRow(${prevIdx}, false)"`}>&larr; Prev Hit</button>
        <button type="button" class="btn-secondary" ${nextIdx < 0 ? 'disabled style="opacity:0.45;"' : `onclick="openBatchRow(${nextIdx}, false)"`}>Next Hit &rarr;</button>
        <button type="button" class="btn-secondary" onclick="scrollToBatchMatrix()">&uarr; Back to Matrix</button>
      </div>
    </div>`;
}

function scrollToBatchMatrix() {
  const el = document.getElementById('batch-results-card');
  if (el) el.scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function openBatchRow(idx, shouldScroll = true) {
  const row = BATCH_ROWS[idx];
  if (!row || !row.payload) return;
  SELECTED_BATCH_IDX = idx;
  const sp = row.payload.report.species_id || getScopeSpecies('batch');
  const tissue = sp === 'human' ? document.getElementById('batch-tissue').value.trim() : '';
  const cell = document.getElementById('batch-cell').value.trim();
  syncUrlPermalink(row.payload.report.gene_symbol, tissue, cell, sp);
  renderBatchTable();
  applyLoadedPayload(row.payload);
  renderBatchReportNav();
  if (shouldScroll) {
    document.getElementById('results-area').scrollIntoView({ behavior: 'smooth', block: 'start' });
  }
}

function exportBatchCsv() {
  const ready = BATCH_ROWS.filter(r => r.payload);
  if (!ready.length) return;
  const lines = ['Species,Assembly,Gene,Gene_ID,Risk_Level,Recommended_Modality,Chronos,CCLE_Genetics,Tissue_TPM,Top_Paralog,Best_Exon,Top_Spacer_20nt,Top_U6_Fwd_Oligo,Top_U6_Rev_Oligo'];
  ready.forEach(({ payload }) => {
    const r = payload.report;
    const ess = r.essentiality;
    const chronos = (ess.matched_cell_line && ess.matched_cell_line_chronos !== null) ? ess.matched_cell_line_chronos : (ess.depmap_median_chronos ?? '');
    const clg = r.cell_line_genetics;
    const ccle = (clg && clg.matched_ccle_sample) ? `${(clg.mutations || []).join(';') || 'WT'} | ${clg.cna_label || ''}` : '';
    const tpm = r.expression.requested_tissue_tpm ?? '';
    const par = (r.paralogs && r.paralogs[0]) ? `${r.paralogs[0].gene_symbol} (${r.paralogs[0].buffering_risk})` : 'None';
    const bestEx = r.exons.filter(e => e.is_coding).sort((a, b) => b.priority_score - a.priority_score)[0];
    const topG = (r.recommended_guides && r.recommended_guides[0]) || (r.crispri_guides && r.crispri_guides[0]);
    lines.push([
      r.species_common_name || 'Human',
      r.assembly,
      r.gene_symbol,
      r.ensembl_id,
      r.modality.risk_level,
      r.modality.primary_modality,
      chronos,
      ccle,
      tpm,
      par,
      bestEx ? `Exon${bestEx.exon_number} (${bestEx.priority_score}/100)` : '',
      topG ? topG.spacer_20nt : '',
      topG ? topG.cloning_fwd_oligo : '',
      topG ? topG.cloning_rev_oligo : ''
    ].map(csvCell).join(','));
  });
  triggerCsvDownload('crispr_preflight_batch_triage.csv', lines);
}

window.addEventListener('DOMContentLoaded', () => {
  renderExampleChips();
  const urlParams = new URLSearchParams(window.location.search);
  const qSpecies = (urlParams.get('species') || '').trim().toLowerCase();
  if (qSpecies && SPECIES_META[qSpecies]) {
    document.getElementById('species-input').value = qSpecies;
    document.getElementById('batch-species').value = qSpecies;
    onSpeciesChange('single', false);
    onSpeciesChange('batch', false);
  }
  const qGene = urlParams.get('gene');
  const qTissue = urlParams.get('tissue');
  const qCell = urlParams.get('cell_line');
  if (qTissue !== null && getScopeSpecies('single') === 'human') {
    document.getElementById('tissue-input').value = qTissue;
  }
  if (qCell !== null) setCellValue('single', qCell);
  if (qGene && qGene.trim()) {
    document.getElementById('gene-input').value = qGene.trim();
    runCheck(null, false);
  }
});
</script>
</body>
</html>
"""


def make_handler(use_cache: bool = True) -> type[BaseHTTPRequestHandler]:
    """Create an HTTP request handler bound to a shared thread-safe PreflightRunner."""
    shared_http = CachedHttpClient(use_cache=use_cache)
    shared_runner = PreflightRunner(http=shared_http, check_pams=True)
    hot_lock = threading.Lock()
    hot_payloads: dict[tuple[str, str, Optional[str], Optional[str]], tuple[bytes, float]] = {}

    class PreflightRequestHandler(BaseHTTPRequestHandler):
        def log_message(self, fmt: str, *args: Any) -> None:
            # Keep server console logs clean and readable
            return

        def _send_raw_json(self, status: int, raw: bytes) -> None:
            self.send_response(status)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Cache-Control", "no-store, max-age=0")
            self.send_header("Content-Length", str(len(raw)))
            self.end_headers()
            self.wfile.write(raw)

        def _send_json(self, status: int, payload: dict[str, Any]) -> None:
            self._send_raw_json(status, json.dumps(payload).encode("utf-8"))

        def do_POST(self) -> None:  # noqa: N802
            parsed = urllib.parse.urlparse(self.path)
            qs = urllib.parse.parse_qs(parsed.query)
            route_hint = (qs.get("_route") or qs.get("__route") or [""])[0].strip()
            post_path = f"/api/{route_hint}" if route_hint in ("snp-check",) else parsed.path
            if post_path == "/api/snp-check":
                try:
                    length = int(self.headers.get("Content-Length") or 0)
                except ValueError:
                    length = -1
                if length <= 0 or length > 65536:
                    self._send_json(HTTPStatus.BAD_REQUEST, {"error": "Invalid request payload size"})
                    return
                try:
                    raw_body = self.rfile.read(length)
                    payload = json.loads(raw_body.decode("utf-8"))
                except Exception:
                    self._send_json(HTTPStatus.BAD_REQUEST, {"error": "Malformed JSON request body"})
                    return
                specs = payload.get("specs") if isinstance(payload, dict) else None
                if not isinstance(specs, list):
                    self._send_json(HTTPStatus.BAD_REQUEST, {"error": "Missing 'specs' list"})
                    return
                try:
                    results = shared_runner.cellline_variants.check_guides_common_snps(specs)
                except Exception:
                    results = {}
                self._send_json(HTTPStatus.OK, {"results": results})
                return

            self._send_json(HTTPStatus.NOT_FOUND, {"error": "Not found"})

        def do_GET(self) -> None:  # noqa: N802
            parsed = urllib.parse.urlparse(self.path)
            qs = urllib.parse.parse_qs(parsed.query)
            route_hint = (qs.get("_route") or qs.get("__route") or [""])[0].strip()
            path = f"/api/{route_hint}" if route_hint in ("check", "snp-check", "health") else parsed.path

            if path in ("/", "/index.html"):
                body = WEB_APP_HTML.encode("utf-8")
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Cache-Control", "no-store, no-cache, must-revalidate, max-age=0")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return

            if path == "/api/health":
                try:
                    entries = int(shared_http.cache_stats().get("entries", 0))
                except Exception:
                    entries = 0
                self._send_json(
                    HTTPStatus.OK,
                    {
                        "status": "ok",
                        "version": __version__,
                        "cache": {"enabled": use_cache, "entries": entries},
                    },
                )
                return

            if path == "/api/snp-check":
                qs = urllib.parse.parse_qs(parsed.query)
                raw_specs = (qs.get("specs") or [""])[0]
                try:
                    specs = json.loads(raw_specs) if raw_specs else []
                except Exception:
                    self._send_json(HTTPStatus.BAD_REQUEST, {"error": "Invalid specs JSON"})
                    return
                if not isinstance(specs, list):
                    self._send_json(HTTPStatus.BAD_REQUEST, {"error": "specs must be a JSON array"})
                    return
                try:
                    results = shared_runner.cellline_variants.check_guides_common_snps(specs)
                except Exception:
                    results = {}
                self._send_json(HTTPStatus.OK, {"results": results})
                return

            if path == "/api/check":
                qs = urllib.parse.parse_qs(parsed.query)
                species = (qs.get("species") or ["human"])[0].strip() or "human"
                gene = (qs.get("gene") or [""])[0].strip()
                tissue = (qs.get("tissue") or [""])[0].strip() or None
                cell_line = (qs.get("cell_line") or [""])[0].strip() or None
                if not gene:
                    self._send_json(
                        HTTPStatus.BAD_REQUEST,
                        {"error": "Please provide a gene symbol (e.g. STAT3, Trp53, smarca4a, Stat92E, cep-1)."},
                    )
                    return

                hot_key = (
                    species.lower(),
                    gene.upper(),
                    tissue.lower() if tissue else None,
                    cell_line.lower() if cell_line else None,
                )
                if use_cache:
                    with hot_lock:
                        hit = hot_payloads.get(hot_key)
                    if hit is not None and (time.time() - hit[1]) < 900.0:
                        self._send_raw_json(HTTPStatus.OK, hit[0])
                        return

                try:
                    report = shared_runner.run(
                        gene_query=gene,
                        tissue=tissue,
                        cell_line=cell_line,
                        species=species,
                    )
                except ValueError as exc:
                    self._send_json(HTTPStatus.BAD_REQUEST, {"error": str(exc)})
                    return
                except NetworkFetchError:
                    self._send_json(
                        HTTPStatus.BAD_GATEWAY,
                        {
                            "error": (
                                "An upstream public biological database timed out or is temporarily "
                                "unavailable. Please try again in a few seconds."
                            )
                        },
                    )
                    return
                except Exception:
                    self._send_json(
                        HTTPStatus.INTERNAL_SERVER_ERROR,
                        {"error": "Unexpected error while running biological flight check."},
                    )
                    return

                raw = json.dumps(
                    {
                        "report": report.to_dict(),
                        "standalone_html": render_html(report),
                        "markdown": render_markdown(report),
                    }
                ).encode("utf-8")
                if use_cache:
                    with hot_lock:
                        if len(hot_payloads) >= 256:
                            hot_payloads.clear()
                        hot_payloads[hot_key] = (raw, time.time())
                self._send_raw_json(HTTPStatus.OK, raw)
                return

            self._send_json(HTTPStatus.NOT_FOUND, {"error": "Not found"})

    return PreflightRequestHandler


def check_live_dns() -> None:
    """Verify outbound DNS resolution works so live uncached gene lookups succeed."""
    socket.getaddrinfo("www.ebi.ac.uk", 443)


def run_web_server(
    host: str = "127.0.0.1",
    port: int = 7860,
    use_cache: bool = True,
    verify_network: bool = False,
) -> None:
    """Start the multi-threaded crispr-preflight web server."""
    if verify_network:
        check_live_dns()
    handler_cls = make_handler(use_cache=use_cache)
    server = ThreadingHTTPServer((host, port), handler_cls)
    print(f"[*] crispr-preflight web app running at http://{host}:{port}", flush=True)
    print("    Press Ctrl+C to stop.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()

"use strict";
// BOLDTrace web page. DOM is built with textContent only (no innerHTML).

const state = { options: null, token: null, job: null, assignments: [], references: [] };

function el(tag, attrs = {}, text) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") e.className = v; else e.setAttribute(k, v);
  }
  if (text !== undefined && text !== null) e.textContent = String(text);
  return e;
}
const $ = (id) => document.getElementById(id);

async function api(path, body) {
  const opts = body === undefined ? {} :
    { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
  const res = await fetch(path, opts);
  const data = await res.json().catch(() => ({ error: "invalid response" }));
  if (!res.ok) throw new Error(data.error || ("HTTP " + res.status));
  return data;
}

function fmt(v, digits) {
  if (v === null || v === undefined || v === "") return "";
  if (typeof v === "number") return digits === undefined ? String(v) : v.toFixed(digits);
  if (typeof v === "boolean") return v ? "yes" : "no";
  return String(v);
}

function kpis(container, items) {
  container.replaceChildren();
  for (const [k, v] of items) {
    const box = el("div");
    box.append(el("span", { class: "kpi-value" }, fmt(v)), el("span", { class: "kpi-label" }, k));
    container.append(box);
  }
}

// Sortable table. columns: [{key, label, num, digits}]
function renderTable(table, rows, columns, onClick, rowClass) {
  table.replaceChildren();
  const thead = el("thead"), tr = el("tr");
  let sortKey = null, asc = true;
  for (const c of columns) {
    const th = el("th", { class: c.num ? "num" : "" }, c.label);
    th.addEventListener("click", () => {
      asc = sortKey === c.key ? !asc : true; sortKey = c.key;
      const sorted = [...rows].sort((a, b) => {
        const x = a[c.key], y = b[c.key];
        if (x === y) return 0;
        if (x === null || x === undefined) return 1;
        if (y === null || y === undefined) return -1;
        return (x < y ? -1 : 1) * (asc ? 1 : -1);
      });
      body(sorted);
    });
    tr.append(th);
  }
  thead.append(tr);
  const tbody = el("tbody");
  function body(list) {
    tbody.replaceChildren();
    for (const r of list) {
      const row = el("tr", { class: (onClick ? "clickable " : "") + (rowClass ? rowClass(r) : "") });
      for (const c of columns) row.append(el("td", { class: c.num ? "num" : "" }, fmt(r[c.key], c.digits)));
      if (onClick) row.addEventListener("click", () => {
        tbody.querySelectorAll("tr.selected").forEach((x) => x.classList.remove("selected"));
        row.classList.add("selected"); onClick(r);
      });
      tbody.append(row);
    }
  }
  body(rows);
  table.append(thead, tbody);
}

function showTab(name) {
  document.querySelectorAll("#tabs button").forEach((b) =>
    b.setAttribute("aria-selected", b.dataset.tab === name ? "true" : "false"));
  document.querySelectorAll(".tab").forEach((s) => { s.hidden = s.id !== "tab-" + name; });
  if (name === "stats") loadStats();
}

function fillList(ul, items) {
  ul.replaceChildren();
  for (const i of items) ul.append(el("li", {}, i));
}

async function loadOptions() {
  const o = await api("api/options");
  state.options = o;
  const dbSel = $("db-select"), modeSel = $("mode-select"), acc = $("accept-select");
  for (const d of o.databases) {
    const opt = el("option", { value: d.key }, d.key.replace(/_/g, " ") + " (" + d.code + ")");
    if (d.key === "animal_public") opt.selected = true;
    dbSel.append(opt);
  }
  for (const m of o.modes) {
    modeSel.append(el("option", { value: m.key },
      m.key.replace(/_/g, " ") + " — min " + Math.round(m.mi * 100) + "%, ≤" + m.maxh + " hits"));
  }
  for (const r of o.ranks) {
    const opt = el("option", { value: r }, r);
    if (r === o.accept_rank) opt.selected = true;
    acc.append(opt);
  }
  const f = $("analysis-form");
  for (const k of ["min_identity", "min_coverage", "max_delta", "max_matches", "min_align_length", "min_agreement"]) {
    f.elements[k].value = o.defaults[k];
  }
  f.elements.rank_cap.checked = o.defaults.rank_cap;
  const grid = $("threshold-grid");
  for (const r of o.ranks) {
    const lab = el("label", {}, r);
    const inp = el("input", { class: "text-input", type: "number", min: "0", max: "100", step: "0.1",
                              name: "thr_" + r, "aria-label": "Minimum identity for " + r });
    inp.value = o.defaults.rank_thresholds[r];
    lab.append(inp); grid.append(lab);
  }
  const syncCap = () => { $("thresholds").disabled = !f.elements.rank_cap.checked; };
  f.elements.rank_cap.addEventListener("change", syncCap); syncCap();
  $("global-note").textContent = o.global_note;
}

async function upload() {
  const file = $("fasta-file").files[0];
  const text = file ? await file.text() : $("fasta-text").value;
  const filename = file ? file.name : "pasted.fasta";
  $("upload-btn").disabled = true;
  try {
    const r = await api("api/upload", { fasta: text, filename });
    state.token = r.token;
    state.job = null; state.jobInfo = null;
    $("reanalysis-box").hidden = true;
    $("stat-input").textContent = filename;
    $("stat-unique").textContent = r.n_unique;
    $("stat-accepted").textContent = "—"; $("stat-refs").textContent = "—";
    $("upload-result").hidden = false;
    kpis($("upload-kpis"), [["sequences", r.n_records], ["valid", r.n_valid],
      ["unique sequences", r.n_unique], ["duplicates", r.n_duplicates],
      ["errors", r.n_errors], ["warnings", r.warnings.length]]);
    fillList($("upload-errors"), r.errors);
    fillList($("upload-warnings"), r.warnings);
    $("upload-errors-box").hidden = r.errors.length === 0;
    $("upload-errors-box").open = r.errors.length > 0;
    $("upload-warnings-box").hidden = r.warnings.length === 0;
    $("upload-hash").textContent = "File SHA-256 " + r.file_sha256 + " · sequence set hash " + r.set_hash;
    $("analyze-btn").disabled = r.n_unique === 0;
    $("analyze-hint").textContent = r.n_unique + " unique sequence(s) ready.";
  } catch (e) {
    $("upload-result").hidden = false;
    kpis($("upload-kpis"), []);
    fillList($("upload-errors"), [e.message]);
    $("upload-errors-box").hidden = false; $("upload-errors-box").open = true;
    $("upload-warnings-box").hidden = true;
    $("analyze-btn").disabled = true;
  } finally {
    $("upload-btn").disabled = false;
  }
}

function formParams() {
  const f = $("analysis-form").elements;
  const num = (k) => Number(f[k].value);
  return { min_identity: num("min_identity"), min_coverage: num("min_coverage"),
           max_delta: num("max_delta"), max_matches: num("max_matches"),
           min_align_length: num("min_align_length"), min_agreement: num("min_agreement"),
           rank_cap: f.rank_cap.checked,
           rank_thresholds: Object.fromEntries(state.options.ranks.map((r) => [r, num("thr_" + r)])) };
}

async function analyze() {
  const f = $("analysis-form").elements;
  const body = {
    token: state.token, project: f.project.value, user: f.user.value,
    db: f.db.value, mode: f.mode.value, accept_rank: f.accept_rank.value,
    params: formParams(),
  };
  $("reanalysis-box").hidden = true;
  $("analyze-btn").disabled = true;
  $("analyze-btn").classList.add("busy");
  $("job-box").hidden = false;
  $("job-error").hidden = true;
  $("job-kpis").replaceChildren();
  try {
    const r = await api("api/analyze", body);
    state.job = r.id;
    poll();
  } catch (e) {
    $("job-status").textContent = "not started";
    $("job-error").textContent = e.message; $("job-error").hidden = false;
    $("analyze-btn").disabled = false;
    $("analyze-btn").classList.remove("busy");
  }
}

async function poll() {
  let j;
  try { j = await api("api/job?id=" + encodeURIComponent(state.job)); }
  catch (e) { $("job-error").textContent = e.message; $("job-error").hidden = false; return; }
  $("job-status").textContent = j.status;
  fillList($("job-messages"), j.messages);
  if (j.status === "queued" || j.status === "running") { setTimeout(poll, 3000); return; }
  $("analyze-btn").disabled = false;
  $("analyze-btn").classList.remove("busy");
  if (j.status === "failed" || j.status === "cancelled") {
    $("job-error").textContent = j.error || "The analysis was cancelled because the server was stopped.";
    $("job-error").hidden = false; return;
  }
  await showFinished(j);
}

// Shows a finished analysis (from BOLD or recomputed) and enables the
// recomputation on its data.
async function showFinished(j) {
  state.job = j.id;
  state.jobInfo = j;
  const s = j.summary;
  kpis($("job-kpis"), [["unique sequences", s.n_unique], ["accepted assignments", s.n_accepted],
    ["assigned to species", s.by_rank.species || 0], ["no admissible match", s.status.no_admissible_match || 0],
    ["no BOLD match", s.status.no_matches || 0], ["references in LCA", s.n_references_in_lca],
    ["process IDs retrieved", s.n_process_ids]]);
  $("stats-project").value = j.project;
  $("stat-accepted").textContent = s.n_accepted;
  $("stat-refs").textContent = s.n_references_in_lca;
  $("reanalysis-box").hidden = false;
  syncReanalysisHint();
  await loadResults();
}

function syncReanalysisHint() {
  const j = state.jobInfo, f = $("analysis-form").elements;
  if (!j) return;
  const changed = f.db.value !== j.db || f.mode.value !== j.mode;
  $("reanalysis-hint").textContent = changed
    ? "The BOLD library or search mode differs from the downloaded data (" + j.db + ", " + j.mode +
      "): recomputing keeps the downloaded data; use Run analysis to query BOLD again."
    : "Applies the parameters above to the BOLD hits and metadata already retrieved (" + j.db + ", " +
      j.mode + "): no new BOLD request, results in seconds.";
}

async function reanalyze() {
  const f = $("analysis-form").elements;
  const btn = $("reanalyze-btn");
  btn.disabled = true; btn.classList.add("busy");
  $("job-error").hidden = true;
  try {
    const j = await api("api/reanalyze", { id: state.job, accept_rank: f.accept_rank.value,
      params: formParams(), record: $("reanalysis-record").checked });
    $("job-status").textContent = j.status;
    fillList($("job-messages"), j.messages);
    await showFinished(j);
  } catch (e) {
    $("job-error").textContent = e.message; $("job-error").hidden = false;
  } finally {
    btn.disabled = false; btn.classList.remove("busy");
  }
}

async function loadResults() {
  const id = encodeURIComponent(state.job);
  const [assignments, references, attribution] = await Promise.all([
    api("api/assignments?id=" + id), api("api/references?id=" + id), api("api/attribution?id=" + id)]);
  state.assignments = assignments; state.references = references;
  document.querySelectorAll("#tabs button:disabled").forEach((b) => { b.disabled = false; });
  drawAssignments();
  drawReferences();
  renderTable($("attr-table"), attribution, [
    { key: "identifier", label: "Taxonomic identifier" }, { key: "orcid", label: "ORCID" },
    { key: "n_records", label: "BOLD records", num: true },
    { key: "n_asvs_supported", label: "ASVs supported", num: true },
    { key: "n_lca_participations", label: "LCA participations", num: true },
    { key: "n_taxa", label: "Distinct taxa", num: true },
    { key: "fractional_credit", label: "Fractional credit (experimental)", num: true, digits: 3 }]);
  const dl = (file) => "api/export?id=" + id + "&file=" + file;
  $("dl-assignments").href = dl("assignments.csv");
  $("dl-matches").href = dl("bold_matches.csv");
  $("dl-references").href = dl("reference_usage.csv");
  $("dl-attribution").href = dl("taxonomic_attribution.csv");
  const list = $("export-list");
  list.replaceChildren();
  for (const name of ["assignments.csv", "bold_matches.csv", "reference_usage.csv",
                      "taxonomic_attribution.csv", "analysis_provenance.json", "used_process_ids.tsv"]) {
    const li = el("li"), a = el("a", { href: "api/export?id=" + id + "&file=" + name, download: name }, name);
    li.append(a); list.append(li);
  }
}

function drawAssignments() {
  $("match-panel").hidden = true;
  const q = $("assign-filter").value.trim().toLowerCase();
  const rows = state.assignments.filter((r) => !q ||
    String(r.asv_id).toLowerCase().includes(q) || String(r.taxon || "").toLowerCase().includes(q));
  renderTable($("assign-table"), rows, [
    { key: "asv_id", label: "ASV ID" }, { key: "taxon", label: "Taxon" }, { key: "rank", label: "Rank" },
    { key: "best_identity", label: "Best identity (%)", num: true, digits: 2 },
    { key: "best_coverage", label: "Coverage (%)", num: true, digits: 0 },
    { key: "n_selected", label: "Matches in LCA", num: true },
    { key: "agreement", label: "Agreement (%)", num: true, digits: 1 },
    { key: "n_dissenting", label: "Dissenting", num: true },
    { key: "n_matches", label: "Matches retrieved", num: true },
    { key: "status", label: "Status" }, { key: "flags", label: "Notes" }], showMatches);
}

async function showMatches(r) {
  const panel = $("match-panel");
  panel.hidden = false;
  const link = $("dl-asv-matches");
  if (!r.sequence_hash) {
    $("match-title").textContent = r.asv_id + ": " + (r.error || "no sequence");
    $("match-table").replaceChildren(); link.hidden = true; return;
  }
  link.hidden = false;
  link.href = "api/export?id=" + encodeURIComponent(state.job) + "&file=bold_matches.csv&hash=" +
    encodeURIComponent(r.sequence_hash);
  $("match-title").textContent = "References for " + r.asv_id + (r.taxon ? " → " + r.taxon + " (" + r.rank + ")" : "");
  const rows = await api("api/matches?id=" + encodeURIComponent(state.job) + "&hash=" + encodeURIComponent(r.sequence_hash));
  renderTable($("match-table"), rows, [
    { key: "process_id", label: "Process ID" }, { key: "marker", label: "Marker" },
    { key: "reference_taxon", label: "Taxon" }, { key: "identified_by", label: "Identified by" },
    { key: "identity", label: "Identity (%)", num: true, digits: 2 },
    { key: "coverage", label: "Coverage (%)", num: true, digits: 0 },
    { key: "align_length", label: "Length", num: true },
    { key: "included_in_lca", label: "In LCA" }, { key: "decision", label: "Decision" },
    { key: "dissent_rank", label: "Dissents at" }, { key: "supports_assignment", label: "Supports" },
    { key: "bin_uri", label: "BIN" }], null, (m) => m.supports_assignment ? "" : "excluded");
}

function drawReferences() {
  const q = $("ref-filter").value.trim().toLowerCase();
  const rows = state.references.filter((r) => !q || [r.process_id, r.taxon, r.identified_by]
    .some((v) => String(v || "").toLowerCase().includes(q)));
  renderTable($("ref-table"), rows, [
    { key: "process_id", label: "Process ID" }, { key: "marker", label: "Marker" },
    { key: "taxon", label: "Taxon" }, { key: "identified_by", label: "Identifier" },
    { key: "retrieved_usage", label: "Retrieved", num: true },
    { key: "lca_participation", label: "LCA", num: true },
    { key: "assignment_contribution", label: "Contribution", num: true },
    { key: "recordsets", label: "Projects/datasets" }]);
}

function kv(table, obj, labels) {
  table.replaceChildren();
  for (const [k, label] of labels) {
    const tr = el("tr");
    tr.append(el("th", {}, label), el("td", { class: "num" }, fmt(obj[k])));
    table.append(tr);
  }
}

const STAT_LABELS = [["n_analyses", "Analyses"], ["n_projects", "Distinct projects"],
  ["n_events", "Identification events (ASV × reference)"], ["n_lca_events", "Events in LCA sets"],
  ["n_contribution_events", "Events in accepted assignments"], ["n_asvs", "Distinct ASVs"],
  ["n_process_ids", "Process IDs retrieved"], ["n_process_ids_in_lca", "Process IDs in LCA sets"],
  ["n_reference_sequences", "Distinct reference sequences"]];

async function loadStats() {
  const p = encodeURIComponent($("stats-project").value.trim());
  const u = encodeURIComponent($("stats-user").value.trim());
  try {
    const s = await api("api/stats?project=" + p + "&user=" + u);
    kv($("stats-personal"), s.personal, STAT_LABELS);
    kv($("stats-global"), s.global_, STAT_LABELS);
  } catch (e) { $("global-note").textContent = e.message; }
}

document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("#tabs button").forEach((b) => b.addEventListener("click", () => showTab(b.dataset.tab)));
  $("upload-btn").addEventListener("click", upload);
  $("fasta-file").addEventListener("change", () => {
    const f = $("fasta-file").files[0];
    $("file-name").textContent = f ? f.name : "";
    if (f) $("fasta-text").value = "";
  });
  $("fasta-text").addEventListener("input", () => {
    if ($("fasta-text").value) { $("fasta-file").value = ""; $("file-name").textContent = ""; }
  });
  $("analyze-btn").addEventListener("click", analyze);
  $("reanalyze-btn").addEventListener("click", reanalyze);
  $("db-select").addEventListener("change", syncReanalysisHint);
  $("mode-select").addEventListener("change", syncReanalysisHint);
  $("assign-filter").addEventListener("input", drawAssignments);
  $("ref-filter").addEventListener("input", drawReferences);
  $("stats-btn").addEventListener("click", loadStats);
  loadOptions().catch((e) => { $("analyze-hint").textContent = "Cannot load options: " + e.message; });
});

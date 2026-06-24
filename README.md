# MarkerAtlas — Single-Cell Marker Gene Reference

A browser-based reference tool for single-cell RNA sequencing researchers. Eight modes covering the full annotation workflow — from browsing cell types through batch-annotating all clusters and generating publication-ready visualisation code.

**Live tool:** `https://chakitarora.github.io/markeratlas`

Requires HTTP server (`python3 -m http.server 8080` locally, or GitHub Pages). Does not work as `file://`. No account. No data leaves your machine.

---

## What's new in this version

- **Batch annotation** — annotate all clusters at once from a tab-separated table
- **Alias resolution** — CD45, PD-1, CD56, SMA, TIM-3 and 116 other aliases resolve automatically everywhere
- **Evidence tiers** — Gold/Silver/Bronze badges on every marker (literature support depth)
- **Confidence threshold** — slider to require minimum canonical matches before showing a result
- **Session history** — all cluster annotations logged during the session, click to reload
- **PMID hover cards** — hover any reference link to see title, journal, authors without leaving the tool
- **Threshold filter** — require minimum canonical matches to reduce noise in long gene lists

---

## Repository contents

```
markeratlas/
├── index.html                  Complete tool — 1784 lines, 100 KB
├── marker_db.json              Full database as JSON (83 KB)
├── marker_associations.tsv     339 gene–celltype associations (now includes tier)
├── celltypes.tsv               36 cell type definitions (now includes cl_id, minimal_panel)
└── README.md                   This file
```

---

## Run locally

```bash
cd markeratlas/
python3 -m http.server 8080
# open http://localhost:8080
```

---

## The eight modes

### 1 — Browse

Displays all 36 cell types as cards with canonical (blue), supporting (grey), and negative (rust) marker tags. Filter by **Category** and **Tissue** in the sidebar. Click any gene tag to instantly switch to Gene mode. Click any card to open the full detail panel.

Each gene tag now shows an evidence tier badge:
- **G** (Gold, amber) — validated in >20 independent studies. CD3E, CD68, PECAM1, MBP, ALB, etc.
- **S** (Silver, slate) — validated in 5–20 studies. Most supporting markers.
- **B** (Bronze, brown) — fewer studies or highly context-specific. SPINK2, TREM2, FAP, etc.

---

### 2 — Gene

Type any gene symbol — or an alias. **116 aliases are resolved automatically:**

| Alias | Resolves to | | Alias | Resolves to |
|---|---|---|---|---|
| CD45 | PTPRC | | PD-1 | PDCD1 |
| CD56 | NCAM1 | | PD-L1 | CD274 |
| CD31 | PECAM1 | | TIM-3 | HAVCR2 |
| SMA | ACTA2 | | LAG-3 | LAG3 |
| NeuN | RBFOX3 | | CD25 | IL2RA |
| CD20 | MS4A1 | | CD127 | IL7R |
| VE-Cadherin | CDH5 | | Blimp1 | PRDM1 |
| NG2 | CSPG4 | | IBA1 | AIF1 |

When an alias is resolved, an ochre note shows the mapping (e.g. `CD45 → PTPRC`). If no exact or alias match is found, "Did you mean:" suggestions appear from the alias table.

**Example — typing CD45:**
```
→ Resolved alias: CD45 → PTPRC
PTPRC — 0 cell type associations in MarkerAtlas database
(PTPRC/CD45 is a pan-leukocyte marker not individually indexed per subtype)
```

**Example — typing PD-1:**
```
→ Resolved alias: PD-1 → PDCD1
PDCD1 — 1 cell type association
Exhausted T cell  [canonical]  Tumor
```

---

### 3 — Cluster

Paste top DE genes from one cluster. Now includes:

**Tissue-aware scoring** — select tissue to boost cell types validated there and penalise atypical ones.

**Confidence threshold** — slider (0–4) sets the minimum number of canonical marker matches required. At 0: show all results. At 2: only cell types with ≥2 canonical markers matched. Reduces noise for long gene lists.

**Alias resolution** — aliases in your gene list are resolved automatically. A note shows what was resolved.

**Negative marker check** — rust panel below each result shows which negative markers to verify are absent.

**Session history** — each annotation is logged in the sidebar. Click any entry to reload that gene list.

**Export** — TSV table with rank, cell type, score, matched genes, tissue prior, PMIDs.

---

### 4 — Batch

**Annotate all clusters at once.** Paste a table — one cluster per line — and get a full annotation summary in seconds.

**Input format:**

Tab-separated (cluster ID + gene list):
```
Cluster0    CD3E,CD8A,GZMB,PRF1,NKG7,GNLY
Cluster1    C1QA,C1QB,CD68,MRC1,APOE,TREM2
Cluster2    GFAP,AQP4,S100B,ALDH1L1,VIM
Cluster3    COL1A1,COL1A2,DCN,LUM,PDGFRA
```

Gene-list only (clusters numbered automatically):
```
CD3E,CD8A,GZMB,PRF1,NKG7
C1QA,C1QB,CD68,MRC1,APOE
GFAP,AQP4,S100B,ALDH1L1
```

**Output table columns:** Cluster | Top annotation | Score | Confidence | Matched genes | Alt 1 | Alt 2

**Confidence rating:**
- **High** — top score >2× the second-best candidate
- **Medium** — top score 1.3–2× second-best
- **Low** — top score <1.3× second-best (annotation is ambiguous)

Click any annotated cell type to open its detail panel.

**Export** — full TSV with all columns including alternative candidates, unmatched genes, and tissue prior.

**Seurat integration:**
```r
# Export your FindAllMarkers results, reformat as batch input
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
top20 <- markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
# Paste: cluster<TAB>gene1,gene2,...  into Batch mode
```

```python
# Scanpy
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
result = adata.uns['rank_genes_groups']
for cluster in adata.obs['leiden'].unique():
    genes = result['names'][cluster][:20]
    print(f"{cluster}\t{','.join(genes)}")
```

---

### 5 — Co-exp

Enter 2–6 gene symbols. Returns cell types where **all** are present (canonical or supporting), plus partial matches (≥50%). Aliases resolved automatically.

Each matched gene shows a badge: **C** = canonical, **S** = supporting. Missing genes shown with strikethrough in partial matches.

**Example — CD3E, CD8A, GZMB:**
```
All 3 present (1 cell type):
  CD8+ T cell    CD3E[C]  CD8A[C]  GZMB[C]   Blood · Tumor · Lymph node

Partial (≥2/3):
  NK cell        CD3E[C]  ~~CD8A~~  GZMB[C]   2/3
```

**Example — PD-1, TIM-3 (alias input):**
```
→ PD-1 resolved to PDCD1, TIM-3 resolved to HAVCR2
All 2 present (1 cell type):
  Exhausted T cell    PDCD1[C]  HAVCR2[C]    Tumor
```

Ten pre-loaded example combinations as chips.

---

### 6 — Panel

Select any cell type → get the minimal discriminating gene set, key negative markers, and disambiguation from the most easily confused neighbours.

**Output:**
- **Unique markers** — canonical only in this cell type (most reliable)
- **Semi-unique** — canonical in ≤1 other type
- **Negative markers** — must be absent
- **Disambiguation table** — vs top 3 confused types: shared genes, positive discriminators (green), negative discriminators (rust)

Export per cell type as `.txt`.

---

### 7 — Overlap

Ranked table of all cell type pairs by number of shared markers. Click any cell type name → detail panel. Click any shared gene → gene lookup. Instantly shows which annotation calls are inherently ambiguous.

---

### 8 — DotPlot

After a cluster annotation (mode 3), generates ready-to-run code for:

**Seurat (R):** `RenameIdents()` + `DotPlot()` + optional `FeaturePlot()`
**Scanpy (Python):** `cluster_map`, `marker_dict`, `sc.pl.dotplot()` + UMAP

Select language, set markers per type (2–5), click Generate. Copy button included.

---

## Detail panel

Available from any mode:

| Section | Content |
|---|---|
| Minimal panel | Pre-computed 2–5 gene set |
| Canonical markers | Clickable, runs gene lookup. Tier badge shown. |
| Supporting markers | Clickable |
| Negative markers | Rust — should be absent |
| Notes | Caveats, activation-state changes, co-expression requirements |
| Tissue confidence | Bar chart, 0–100% per tissue |
| Primary references | PubMed links — **hover for title/journal/authors** without leaving the tool |
| Cell Ontology ID | CL: accession, links to OLS4 browser |

---

## Database

**36 cell types · 289 genes · 339 associations · 116 aliases**

Evidence tier per gene:
- **Gold** — CD3E, CD4, CD8A, CD19, MS4A1, PECAM1, CDH5, MBP, SFTPC, ALB, HBB, FOXP3, GFAP, RBFOX3, LGR5, and others validated across hundreds of studies
- **Silver** — most supporting markers and moderately validated canonical markers
- **Bronze** — SPINK2, MLLT3 (HSC-specific), TREM2, SPP1, FAP (context-specific or fewer studies)

New JSON fields:
```json
{
  "gene_index": {
    "CD3E": [{"celltype": "CD4+ T cell", "role": "canonical", "tier": "gold", ...}]
  }
}
```

TSV files updated: `marker_associations.tsv` now includes `tier` column. `celltypes.tsv` now includes `cl_id`, `minimal_panel`, `unique_markers` columns.

---

## Using data files

```r
# R — load with tier
m <- read.delim("marker_associations.tsv")
subset(m, tier == "gold" & role == "canonical")  # gold-tier canonicals only
subset(m, gene == "TREM2")                        # all associations for TREM2
```

```python
# Python
import pandas as pd, json
m = pd.read_csv("marker_associations.tsv", sep="\t")
m[(m.tier == "gold") & (m.role == "canonical")]  # gold canonicals
m[m.gene.isin(["CD3E","CD8A","GZMB"])]          # any of these genes
```

---

## Comparison with related tools

| Feature | MarkerAtlas | ACT | sctype | Azimuth | CellTypist |
|---|---|---|---|---|---|
| Input | Marker list | Marker list | Marker list | Expression matrix | Expression matrix |
| Infrastructure | Browser, offline | Remote server | Browser/R | Remote server | Python |
| Transparency | Full | Score only | Partial | Score only | Score only |
| Alias resolution | Yes (116) | No | No | No | No |
| Tissue-aware scoring | Yes | No | No | Yes | No |
| Negative markers | Yes | No | No | No | No |
| Evidence tiers | Yes | No | No | No | No |
| Batch annotation | Yes | No | No | No | No |
| Co-expression query | Yes | No | No | No | No |
| Minimal panel builder | Yes | No | No | No | No |
| DotPlot code generation | Yes | No | No | No | No |
| PMID hover cards | Yes | No | No | No | No |
| Session history | Yes | No | No | No | No |
| Offline operation | Yes | No | No | No | No |

---

## Limitations

**Species:** Human only.
**Coverage:** 36 common cell types. Rare subtypes not included.
**Cluster annotation:** Scores marker overlap only — not expression levels or trajectory.
**Evidence tiers:** Based on curated estimates, not automated systematic review.
**Dropout:** Negative marker absence may be technical, not biological.
**DotPlot code:** Cluster numbering in generated code corresponds to rank in annotation results — adjust to match your actual cluster IDs.

---

## Disclaimer

Research reference tool only. Not for clinical or diagnostic use. Verify marker suitability for your specific tissue, species, and conditions before use in publications. Cite CellMarker 2.0 and PanglaoDB when using marker information in published work.

---

## Data sources

| Database | Reference |
|---|---|
| CellMarker 2.0 | Zhang X et al., *Nucleic Acids Res.* 2023. [PMID 36300619](https://pubmed.ncbi.nlm.nih.gov/36300619/) |
| PanglaoDB | Franzén O et al., *Database* 2019. [PMID 30951143](https://pubmed.ncbi.nlm.nih.gov/30951143/) |

---

## Deploy

Upload all five files to GitHub, enable Pages. Live at `https://yourusername.github.io/markeratlas`.

## License

MIT

---

*Built by [Chakit Arora](https://chakitarora.github.io) · Post-Doctoral Fellow, BIO@SNS, Scuola Normale Superiore, Pisa*

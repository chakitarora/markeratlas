# MarkerAtlas — Single-Cell Marker Gene Reference

A browser-based reference tool for single-cell RNA sequencing researchers. Eight modes covering the full annotation workflow: browse cell types, look up a gene, annotate one cluster, annotate all clusters at once, query co-expression, build a minimal marker panel, inspect annotation ambiguity, and generate dot plot code for Seurat or Scanpy.

**Live tool:** `https://chakitarora.github.io/markeratlas`

Requires an HTTP server — does not work when opened as a local `file://` URL. Run locally with `python3 -m http.server 8080` then open `http://localhost:8080`. No account. No data leaves your machine.

---

## The gap this fills

Cell type annotation in scRNA-seq is one of the most consequential and least standardised steps in the workflow. The standard approach — run differential expression, Google the top genes, cross-reference CellMarker or PanglaoDB one gene at a time — is fragmented, undocumented, and slow. Researchers typically have 10–20 clusters to annotate, do it manually per cluster, and produce no audit trail for their decisions.

Existing tools fall into two categories. Reference-based tools (Azimuth, CellTypist, SingleR) require either a Python/R environment or upload your expression matrix to a remote server, and return a score without explaining what drove the call. Marker-based lookup tools (CellMarker 2.0, PanglaoDB, ACT, sctype) are designed for one gene or one cluster at a time, do not surface conflicts, and do not help you act on the results once you have them — you still need to write the visualisation code, build the supplementary table, and verify each call by hand.

MarkerAtlas addresses the full workflow: from initial marker lookup through batch annotation of all clusters, to verification with negative markers, to generating the dot plot code for the figure. Everything runs in the browser, everything is transparent, and every association links to its source paper.

---

## Repository contents

```
markeratlas/
├── index.html                  Complete tool — all logic self-contained, 100 KB
├── marker_db.json              Full database as JSON (83 KB, machine-readable)
├── marker_associations.tsv     339 gene–celltype associations as flat TSV
├── celltypes.tsv               36 cell type definitions as flat TSV
└── README.md                   This file
```

---

## The eight modes

---

### Mode 1 — Browse

**What it does:** Displays all 36 cell types as cards. Explore the full marker landscape, understand canonical vs supporting vs negative marker distinctions, and navigate by category or tissue context.

**How to use:**
1. Open the tool — Browse is the default view
2. Filter by **Category** in the left sidebar: Immune, Epithelial, Stromal, Neural, Stem/Progenitor, Tumor microenvironment, Hematopoietic
3. Filter by **Tissue** to show only cell types validated in a specific context: Blood, Brain, Lung, Tumor, Liver, etc.
4. Type in the search box to filter by cell type name, marker gene, or category keyword
5. Click any cell type card to open the full detail panel on the right
6. Click any gene tag on a card to immediately switch to Gene mode for that gene

**What each card shows:**

| Element | Meaning |
|---|---|
| Blue gene tags | Canonical markers — strong, specific, reproducible evidence across multiple studies |
| Grey gene tags | Supporting markers — useful in combination but not definitive alone |
| Rust/red gene tags | Negative markers — should be absent; important for ruling out doublets and misannotation |
| Tissue pills | Top 3 validated tissues with confidence percentage |
| `high` / `medium` badge | Specificity: `high` = reliable across tissues; `medium` = context-dependent |
| G / S / B badge on genes | Evidence tier: Gold (>20 studies), Silver (5–20), Bronze (<5 or highly context-specific) |

**Example — filtering by Neural + Brain:**

Returns 5 cell types. Clicking **Microglia** opens the detail panel:
```
Minimal panel:  CX3CR1  P2RY12  TMEM119
Canonical:      CX3CR1[G]  P2RY12[S]  TMEM119[S]  HEXB[S]  TREM2[B]
Supporting:     SALL1  FCRLS  ITGAM  CD68
Negative:       GFAP  RBFOX3  OLIG2

Notes: P2RY12 and TMEM119 are highly CNS-specific homeostatic markers, lost in
disease-associated microglia (DAM). TREM2 marks the DAM subset — Bronze tier
because it is context-specific and not present in homeostatic microglia.

Tissue confidence:
  Brain   99%  █████████░

Cell Ontology: CL:0000129 ↗
References: PMID 30726726 ↗ (hover for title)  PMID 33558698 ↗
```

---

### Mode 2 — Gene

**What it does:** Takes a single gene symbol — or a common alias — and returns every cell type it marks, with role (canonical or supporting), tissue context, and evidence tier. Essential for interpreting unexpected DE genes and verifying whether a marker is specific or shared across multiple types.

**Alias resolution:** 116 common aliases are resolved automatically. You do not need to know the HGNC-approved symbol:

| You type | Resolves to | You type | Resolves to |
|---|---|---|---|
| CD45 | PTPRC | PD-1 | PDCD1 |
| CD56 | NCAM1 | TIM-3 | HAVCR2 |
| CD20 | MS4A1 | LAG-3 | LAG3 |
| CD31 | PECAM1 | Blimp1 | PRDM1 |
| SMA | ACTA2 | NeuN | RBFOX3 |
| VE-Cadherin | CDH5 | NG2 | CSPG4 |
| CD25 | IL2RA | IBA1 | AIF1 |
| CD127 | IL7R | Tryptase | TPSAB1 |

When resolved, an ochre note confirms the mapping. If nothing is found, "Did you mean:" suggestions appear.

**How to use:**
1. Click **Gene** or click any gene tag in Browse mode
2. Type a gene symbol or alias — results appear immediately
3. Click any cell type name in the results to open its detail panel

Ten example genes are shown as clickable chips, including some aliases to demonstrate the resolution.

**Partial matching:** Typing `SFTPA` returns SFTPB, SFTPC, etc.

**Example 1 — CD45 (alias for pan-leukocyte marker):**
```
→ Resolved alias: CD45 → PTPRC
PTPRC — 0 direct associations
PTPRC is a pan-leukocyte marker not individually indexed per cell type
in the current database. Try CD3E (T cells), CD19 (B cells), CD14 (monocytes).
```

**Example 2 — CD3E (highly specific canonical marker):**
```
CD3E — 3 cell type associations  [Gold tier]

CD4+ T cell       [canonical]   Blood · Lymph node · Tumor
CD8+ T cell       [canonical]   Blood · Tumor · Lymph node
Regulatory T cell [canonical]   Blood · Tumor · Lymph node
```
CD3E alone identifies T cells but cannot distinguish subtypes. Co-expression with CD4, CD8A, or FOXP3 is required.

**Example 3 — ACTA2 (ambiguous multi-type marker):**
```
ACTA2 — 4 cell type associations  [Silver tier]

Smooth muscle cell        [canonical]   Vasculature · Gut · Bladder
Cancer-assoc. fibroblast  [canonical]   Tumor
Fibroblast                [supporting]  Tumor · Skin · Lung
Pericyte                  [supporting]  Brain · Tumor · Kidney
```
ACTA2 (αSMA) marks four distinct cell types. Tissue context is essential — in vasculature it marks smooth muscle, in brain it marks pericytes, in tumor stroma it marks CAFs.

**Example 4 — PD-1 (alias, disease-state marker):**
```
→ Resolved alias: PD-1 → PDCD1
PDCD1 — 1 cell type association  [Silver tier]

Exhausted T cell  [canonical]   Tumor
```
PDCD1 marks disease-associated state. Its expression indicates exhaustion, not a distinct cell type lineage.

**Example 5 — TREM2 (context-specific Bronze marker):**
```
TREM2 — 2 cell type associations  [Bronze tier]

Microglia                [canonical]  Brain
Tumor-assoc. macrophage  [canonical]  Tumor
```
Bronze tier reflects that TREM2 marks disease-associated states in both contexts, not homeostatic identity. Its expression indicates activated/pathological state.

---

### Mode 3 — Cluster

**What it does:** Takes a list of differentially expressed genes from a single cluster and ranks all 36 cell types by marker overlap. Includes tissue-aware scoring, a confidence threshold filter, alias resolution, negative marker verification suggestions, and export.

**How to use:**
1. Click **Cluster**
2. Set **Tissue context** if your sample has a known tissue origin (strongly recommended)
3. Adjust **Min canonical matches** slider if you want to filter out low-evidence results (0 = show all; 2 = require at least 2 canonical markers matched)
4. Paste your top DE genes — one per line or comma-separated. Aliases work.
5. Click **▶ Annotate cluster**

Six pre-loaded example clusters available as chips.

**Where to get input genes:**
```r
# Seurat (R)
markers <- FindMarkers(seurat_obj, ident.1 = 3)
top_genes <- rownames(head(markers[markers$avg_log2FC > 0, ], 20))
cat(top_genes, sep = "\n")
```
```python
# Scanpy (Python)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
top = adata.uns['rank_genes_groups']['names']['3'][:20]
print('\n'.join(top))
```

**Tissue-aware scoring:** Each cell type's raw overlap score is multiplied by its tissue confidence for the selected tissue. Cell types not documented in that tissue are penalised (×0.1). A green badge shows the tissue confidence %; a rust badge flags cell types atypical for your tissue.

**Negative marker check:** Below each ranked candidate, a rust panel shows: "→ Rule out: verify X, Y are absent." These are the negative markers for that cell type not in your input list.

**Session history:** Each annotation is automatically logged in the sidebar history panel. Click any entry to reload that gene list.

**Example — CD8+ T cell cluster, tissue = Blood:**
```
Input: CD3E, CD3D, CD8A, CD8B, GZMB, GZMK, PRF1, NKG7, GNLY, IFNG

1. CD8+ T cell     Blood 99%   10/10 · score 18.0  ████████████
   CD3E CD3D CD8A CD8B GZMB GZMK PRF1 NKG7 GNLY IFNG
   → Rule out: verify CD4, CD19, CD14 are absent

2. NK cell         Blood 98%    4/9  · score 6.9   ████
   GZMB PRF1 NKG7 GNLY
   → Rule out: verify CD3E, CD19, CD14 are absent

3. CD4+ T cell     Blood 99%    2/6  · score 4.0   ██
   CD3E CD3D
```

NK cell at rank 2 is expected — GZMB/PRF1/NKG7/GNLY are shared between cytotoxic T cells and NK cells. CD8A/CD8B are T cell-specific and only appear at rank 1.

**Example — Macrophage cluster from tumor, tissue = Tumor:**
```
Input: C1QA, C1QB, C1QC, CD68, MRC1, APOE, SPP1, TREM2, MARCO, CD163

1. Macrophage               Tumor 95%   8/9 · score 12.4  ██████████
   C1QA C1QB C1QC CD68 MRC1 APOE SPP1 TREM2 MARCO
   → Rule out: verify CD3E, CD19, FCGR3B are absent

2. Tumor-assoc. macrophage  Tumor 95%   3/9 · score 5.7   █████
   SPP1 TREM2 CD68

3. Monocyte                 Tumor —     1/8 · score 0.1   █
   CD68  [penalised: not typical in Tumor]
```

SPP1+TREM2 co-expression at rank 2 suggests a disease-associated or TAM subset. Consider whether this cluster should be split.

**Example — Ambiguous / doublet cluster:**
```
Input: CD3E, COL1A1, GFAP, EPCAM, CD68, PECAM1, HBB, SFTPC

1. Fibroblast      1/9 · score 1
2. Epithelial cell 1/7 · score 1
3. Macrophage      1/9 · score 1
   (all scores equal and very low)
```

No dominant cell type — genes span immune, stromal, epithelial, and neural lineages simultaneously. This strongly suggests doublets. Run DoubletFinder (R) or Scrublet (Python) and re-cluster.

**Export:** TSV table with rank, cell type, score, matched genes, tissue prior, PMIDs — paste directly into a supplementary table.

---

### Mode 4 — Batch

**What it does:** Annotates all clusters at once from a tab-separated table. Produces a full annotation summary with confidence ratings, alternative candidates, and unmatched genes — the core supplementary table for any scRNA-seq paper, generated in one step.

**How to use:**
1. Click **Batch**
2. Set **Tissue context** if known
3. Paste your cluster data — one cluster per line
4. Click **▶ Annotate all clusters**
5. Click **↓ Export full table** to download a TSV

**Input formats:**

Tab-separated (cluster ID + gene list):
```
Cluster0    CD3E,CD8A,GZMB,PRF1,NKG7,GNLY
Cluster1    C1QA,C1QB,CD68,MRC1,APOE,TREM2
Cluster2    GFAP,AQP4,S100B,ALDH1L1,VIM
Cluster3    COL1A1,COL1A2,DCN,LUM,PDGFRA,FAP
Cluster4    SFTPC,SFTPB,ABCA3,ETV5,LAMP3
```

Gene-list only (clusters numbered automatically):
```
CD3E,CD8A,GZMB,PRF1,NKG7
C1QA,C1QB,CD68,MRC1,APOE
GFAP,AQP4,S100B,ALDH1L1
```

Aliases are resolved in gene lists. PD-1, CD45, SMA, TIM-3 etc. all work.

**Output table:** Cluster | Top annotation | Score | Confidence | Matched genes | Alt 1 | Alt 2

**Confidence rating:**
- **High** — top candidate score >2× the second-best
- **Medium** — 1.3–2× second-best
- **Low** — <1.3× second-best; annotation is ambiguous and needs manual review

**Example output:**
```
Cluster0  CD8+ T cell          18.0  High    CD3E, CD8A, GZMB, PRF1, NKG7    NK cell        CD4+ T cell
Cluster1  Macrophage           12.4  High    C1QA, C1QB, CD68, MRC1, APOE    TAM            Monocyte
Cluster2  Astrocyte             9.8  High    GFAP, AQP4, S100B, ALDH1L1      —              —
Cluster3  CAF                   7.2  Medium  COL1A1, FAP, PDGFRA              Fibroblast     —
Cluster4  Alveolar type II      8.0  High    SFTPC, SFTPB, ABCA3, LAMP3      —              —
```

**Seurat workflow:**
```r
# Export all cluster markers
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)
top20 <- all_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 20) %>%
  summarise(genes = paste(gene, collapse = ","))
# Paste top20 as: Cluster0<TAB>gene1,gene2,...
```

```python
# Scanpy workflow
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
for cluster in sorted(adata.obs['leiden'].unique()):
    genes = adata.uns['rank_genes_groups']['names'][cluster][:20]
    print(f"{cluster}\t{','.join(genes)}")
```

Clicking any annotated cell type in the result table opens its detail panel for inspection.

---

### Mode 5 — Co-exp

**What it does:** Given 2–6 gene symbols, returns cell types where **all** of those genes are simultaneously present as canonical or supporting markers. Answers the question: "What cell type co-expresses exactly this combination?" — which neither standard gene lookup nor cluster annotation addresses directly.

**How to use:**
1. Click **Co-exp**
2. Type 2–6 gene symbols, comma or space separated. Aliases work.
3. Press Enter or click **▶ Query**

Ten pre-loaded example combinations as chips.

**Result categories:**
- **Full matches** — all queried genes present. Each gene tagged with **C** (canonical) or **S** (supporting).
- **Partial matches** — ≥50% of genes present. Missing genes shown with strikethrough.

**Example 1 — CD3E, CD8A (T cell co-expression check):**
```
All 2 genes present (2 cell types):
  CD8+ T cell       CD3E[C]  CD8A[C]    Blood · Tumor · Lymph node
  Regulatory T cell CD3E[C]  ~~CD8A~~   — wait, CD8A not in Treg → partial only

Partial (1/2):
  Regulatory T cell  CD3E[C]  ~~CD8A~~
  CD4+ T cell        CD3E[C]  ~~CD8A~~
```

**Example 2 — FOXP3, IL2RA (Treg-specific co-expression):**
```
All 2 genes present (1 cell type):
  Regulatory T cell    FOXP3[C]  IL2RA[C]    Blood · Tumor · Lymph node
```
This combination is uniquely specific to Tregs in the database. Confirmatory for a clean Treg cluster.

**Example 3 — PD-1, TIM-3 (alias input, exhaustion check):**
```
→ PD-1 resolved to PDCD1, TIM-3 resolved to HAVCR2

All 2 genes present (1 cell type):
  Exhausted T cell    PDCD1[C]  HAVCR2[C]    Tumor
```

**Example 4 — SPP1, TREM2 (TAM vs homeostatic macrophage check):**
```
All 2 genes present (2 cell types):
  Tumor-assoc. macrophage    SPP1[C]  TREM2[C]    Tumor
  Macrophage                 SPP1[S]  TREM2[C]    Lung · Liver · Tumor
```
Both cell types co-express these markers. The presence of both in your cluster suggests a disease/tumor-associated macrophage population; tissue context distinguishes them.

**Use cases:**
- Verify a proposed co-expression before committing to an annotation
- Check if a specific gene combination is unique to one cell type or shared
- Identify which cell types could explain unexpected co-expression in your UMAP
- Use with aliases to quickly test clinical marker combinations (CD25/IL2RA + FOXP3 for Tregs)

---

### Mode 6 — Panel

**What it does:** For any cell type, returns the minimal gene set to uniquely identify it, the key negative markers to verify, and a disambiguation guide showing how to distinguish it from the most easily confused neighbours. Useful for flow cytometry panel design, smFISH probe selection, and writing methods sections.

**How to use:**
1. Click **Panel**
2. Select a cell type from the grouped dropdown
3. Click **Build panel**
4. Click **↓ Export panel** to download a `.txt` summary with all markers and PMIDs

**Output elements:**

**Minimal panel** — 2–5 genes labelled `unique` (canonical only in this cell type across the full database) or `semi-unique` (canonical in ≤1 other type). These carry the most discriminating information. Click any gene to run a gene lookup.

**Negative markers** — genes that should be absent. Do not include these in a positive staining panel — use them as exclusion criteria. Note: dropout may mimic true biological absence.

**Disambiguation table** — the 3 most easily confused cell types, with:
- Shared markers (the source of potential confusion)
- Positive discriminators (green +) — present in your type, absent in the confused type
- Negative discriminators (rust −) — absent in your type, present in the confused type

**Example — Regulatory T cell:**
```
Minimal panel:
  FOXP3 [unique]   IL2RA [unique]   IKZF2 [unique]   CTLA4 [semi-unique]

Negative markers:  CD8A   IL7R

Disambiguation:
  vs CD4+ T cell      shared: CD3E, CD4
                       + FOXP3, IL2RA, IKZF2
                       − IL7R (negative in Treg, positive in CD4+ T cell)

  vs Exhausted T cell  shared: TIGIT, CTLA4
                       + FOXP3, IL2RA

  vs CD8+ T cell       shared: CD3E
                       − CD8A (negative in Treg, canonical in CD8+ T)
```

**Example — Oligodendrocyte precursor cell (OPC):**
```
Minimal panel:
  PDGFRA [unique]   CSPG4 [unique]   SOX10 [semi-unique]   OLIG2 [semi-unique]

Negative markers:  MBP   MOG   GFAP   RBFOX3

Disambiguation:
  vs Oligodendrocyte   shared: SOX10, OLIG2, OLIG1
                        + PDGFRA, CSPG4 (OPC-specific)
                        − MBP, MOG (absent in OPC, canonical in mature oligo)
```

PDGFRA and CSPG4 (NG2) are unique to OPC across the full database, making them the most reliable discriminators.

---

### Mode 7 — Overlap

**What it does:** Shows the 20 most ambiguous cell type pairs ranked by number of shared markers. Tells you which annotation calls are inherently harder with marker genes alone — before you start annotating, so you know where to expect Low confidence results and where to look for additional discriminating evidence.

**How to use:** Click **Overlap** — the table renders immediately. No input required. Click any cell type name to open its detail panel. Click any shared gene to run a gene lookup.

**Example output (top rows):**
```
Macrophage          vs  Tumor-assoc. macrophage   C1QA CD68 MRC1 TREM2 APOE ...  8
CD8+ T cell         vs  NK cell                    GZMB PRF1 NKG7 GNLY            4
Fibroblast          vs  Cancer-assoc. fibroblast   THY1 ACTA2 FAP POSTN           4
Epithelial cell     vs  Cholangiocyte              KRT7 KRT19 EPCAM                3
CD4+ T cell         vs  CD8+ T cell                CD3D CD3E                      2
Regulatory T cell   vs  Exhausted T cell           TIGIT CTLA4                    2
Macrophage          vs  Microglia                  TREM2 CD68                      2
Pericyte            vs  Smooth muscle cell         DES ACTA2                       2
```

**Interpretation:**
- Macrophage vs TAM is the hardest call in this database (8 shared markers). A Low confidence result when annotating a macrophage-like cluster is expected — additional markers, tissue context, and functional data are needed to distinguish them.
- CD8+ T cell vs NK cell is the classic cytotoxic ambiguity. CD8A/CD8B (T cell-specific) and NCAM1 (NK-specific) are the key discriminators.
- Knowing these overlaps going in prevents overconfident annotation.

---

### Mode 8 — DotPlot

**What it does:** After running a cluster annotation (mode 3), generates ready-to-run code for a `DotPlot()` in Seurat (R) or `sc.pl.dotplot()` in Scanpy (Python), using the top canonical markers from your annotation. Eliminates the tedious step of assembling the marker list for the standard annotation figure in every scRNA-seq paper.

**How to use:**
1. Run a **Cluster** annotation first (mode 3)
2. Click **DotPlot**
3. Choose language (Seurat / Scanpy) and number of markers per cell type (2–5)
4. Click **Generate code**
5. Click **Copy code** — paste directly into your R or Python session

**What the generated code includes:**

*Seurat output includes:*
- `new.cluster.ids` vector with annotated cell type names
- `RenameIdents()` call to apply the annotation
- `DotPlot()` call with the marker genes, colour scheme, and `RotatedAxis()`
- Optional `FeaturePlot()` for top markers

*Scanpy output includes:*
- `cluster_map` dictionary mapping cluster IDs to cell type names
- `marker_dict` grouped by cell type, ready for `sc.pl.dotplot()`
- `sc.pl.dotplot()` call with `standard_scale="var"` and dendrogram
- UMAP coloured by cell type

**Example — after annotating a CD8+ T cell cluster, Seurat output:**
```r
# Seurat DotPlot — generated by MarkerAtlas
# Cluster annotation: Blood

new.cluster.ids <- c(
  "CD8+ T cell",   # cluster 0
  "NK cell",       # cluster 1
  "CD4+ T cell",   # cluster 2
  ...
)
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

markers_to_plot <- c("CD8A", "CD8B", "GZMB", "NKG7", "PRF1", ...)

DotPlot(
  seurat_obj,
  features = markers_to_plot,
  cols = c("lightgrey", "#2952a3"),
  dot.scale = 8
) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 8))
```

**Note on cluster numbering:** The generated code assigns cell type names in the order of your annotation results (rank 1, 2, 3...). Adjust the `new.cluster.ids` vector to match your actual Seurat or Scanpy cluster IDs before running.

---

## Detail panel

Clicking any cell type name from any mode opens a sliding detail panel on the right:

| Section | Content |
|---|---|
| **Minimal panel** | Pre-computed 2–5 gene set (unique and semi-unique canonical markers) |
| **Canonical markers** | Blue clickable tags — each runs a gene lookup. Evidence tier badge shown. |
| **Supporting markers** | Grey clickable tags — each runs a gene lookup |
| **Negative markers** | Rust tags — should be absent |
| **Notes** | Caveats, subtypes, activation-state changes, co-expression requirements |
| **Tissue confidence** | Bar chart — each validated tissue with 0–100% confidence score |
| **Primary references** | PubMed links — **hover for title, journal, authors** without leaving the tool |
| **Cell Ontology ID** | CL: accession with link to OLS4 browser |

---

## Database

**36 cell types · 289 unique marker genes · 339 gene–celltype associations · 116 aliases**
Human (*Homo sapiens*) only. Sources: CellMarker 2.0 (2023) and PanglaoDB (2022).

### Cell types by category

| Category | Cell types |
|---|---|
| Immune (11) | CD4+ T cell, CD8+ T cell, Regulatory T cell, NK cell, B cell, Plasma cell, Monocyte, Macrophage, Dendritic cell, Mast cell, Neutrophil |
| Epithelial (8) | Epithelial cell, Club cell, Alveolar type I, Alveolar type II, Hepatocyte, Cholangiocyte, Enterocyte, Goblet cell |
| Stromal (5) | Fibroblast, Endothelial cell, Pericyte, Smooth muscle cell, Adipocyte |
| Neural (5) | Neuron, Astrocyte, Oligodendrocyte, Oligodendrocyte precursor cell, Microglia |
| Stem/Progenitor (2) | Intestinal stem cell, Hematopoietic stem cell |
| Tumor microenvironment (3) | Cancer-associated fibroblast, Tumor-associated macrophage, Exhausted T cell |
| Hematopoietic (2) | Erythrocyte, Megakaryocyte |

### Evidence tiers

**Gold** — validated in >20 independent studies, listed in both CellMarker and PanglaoDB, with multiple primary literature PMIDs. Includes: CD3E, CD4, CD8A, CD19, MS4A1, PECAM1, CDH5, MBP, MOG, SFTPC, SFTPB, ALB, HBB, HBA1, FOXP3, GFAP, AQP4, RBFOX3, LGR5, PF4, ADIPOQ, EPCAM, KRT18, KRT19, and others.

**Silver** — validated in 5–20 studies. Includes most supporting markers and moderately validated canonical markers. GZMB, PRF1, NKG7, P2RY12, TMEM119, AGER, PDPN, RGS5, MYH11, ACTA2, CX3CR1, TPSAB1, FABP1, and others.

**Bronze** — fewer studies or highly context-specific. SPINK2, MLLT3 (HSC-specific, few studies), TREM2, SPP1 (disease-state markers, not homeostatic), FAP, IKZF2, CSPG4, PDCD1, HAVCR2.

### JSON schema (marker_db.json)

```json
{
  "celltypes": {
    "CD8+ T cell": {
      "category": "Immune",
      "canonical_markers": ["CD3E", "CD3D", "CD8A", "CD8B", "GZMK", "GZMB"],
      "supporting_markers": ["PRF1", "NKG7", "GNLY", "IFNG"],
      "negative_markers": ["CD4", "CD19", "CD14"],
      "tissue_context": {"Blood": 0.99, "Tumor": 0.95, "Lymph node": 0.97},
      "specificity": "high",
      "notes": "...",
      "pmids": ["30726726", "33558698", "32103181"],
      "sources": ["CellMarker", "PanglaoDB"],
      "cl_id": "CL:0000625",
      "minimal_panel": ["CD8A", "CD8B", "GZMK"],
      "unique_markers": ["CD8A", "CD8B"]
    }
  },
  "gene_index": {
    "CD8A": [{"celltype": "CD8+ T cell", "role": "canonical",
              "confidence": 0.95, "tier": "gold", "tissue_context": {...}}]
  }
}
```

### Using the data files

```r
# R — marker associations with tier
m <- read.delim("marker_associations.tsv")
subset(m, tier == "gold" & role == "canonical")           # highest-confidence markers only
subset(m, gene == "TREM2")                                 # all associations for TREM2
subset(m, category == "Immune" & role == "canonical")      # all immune canonical markers

# Cell type definitions with minimal panels
cts <- read.delim("celltypes.tsv")
subset(cts, category == "Neural")                          # neural cell types
# Parse minimal_panel column (semicolon-separated)
cts$panel_list <- strsplit(cts$minimal_panel, ";")
```

```python
# Python
import pandas as pd, json

m = pd.read_csv("marker_associations.tsv", sep="\t")
m[m.gene == "FOXP3"]                                       # FOXP3 associations
m[(m.tier == "gold") & (m.role == "canonical")]            # gold tier canonicals
m[m.gene.isin(["CD3E","CD8A","GZMB"])]                    # lookup multiple genes

with open("marker_db.json") as f:
    db = json.load(f)
# Minimal panels for all cell types
panels = {ct: data["minimal_panel"] for ct, data in db["celltypes"].items()}
# Cell types with high specificity
high_spec = [ct for ct, data in db["celltypes"].items() if data["specificity"] == "high"]
```

---

## Comparison with related tools

The tools most similar to MarkerAtlas in approach are ACT (xteam.xbio.top/ACT) and sctype (sctype.app), both of which take marker gene lists as input. Azimuth and CellTypist take expression matrices rather than marker lists and are better for automated large-scale annotation. The table below reflects this.

| Feature | MarkerAtlas | ACT | sctype | Azimuth | CellTypist |
|---|---|---|---|---|---|
| **Input type** | Marker gene list | Marker gene list | Marker gene list | Expression matrix | Expression matrix |
| **Infrastructure** | Browser, fully offline | Remote server (uptime-dependent) | Browser / R package | Remote server | Python package |
| **Transparency** | Full — shows which genes matched and why | Score only, no breakdown | Partial | Score only | Score only |
| **Alias resolution** | Yes — 116 aliases (CD45, PD-1, CD56…) | No | No | No | No |
| **Tissue-aware scoring** | Yes — prior multiplier per tissue | No | No | Yes (reference-based) | No |
| **Negative markers** | Yes — shown and checked per result | No | No | No | No |
| **Evidence tiers** | Yes — Gold/Silver/Bronze per gene | No | No | No | No |
| **Batch annotation** | Yes — all clusters in one paste | No | No | No | No |
| **Co-expression query** | Yes — multi-gene intersection | No | No | No | No |
| **Minimal panel builder** | Yes — unique/semi-unique with disambiguation | No | No | No | No |
| **Overlap/confusion matrix** | Yes — most ambiguous pairs ranked | No | No | No | No |
| **DotPlot code generation** | Yes — Seurat + Scanpy | No | No | No | No |
| **PMID hover cards** | Yes — inline paper details | No | No | No | No |
| **Session history** | Yes — all annotations logged | No | No | No | No |
| **Verifiable PMIDs** | Yes — every association linked | No | No | No | No |
| **Cell Ontology IDs** | Yes — CL: per cell type | Partial | No | Yes | Yes |
| **Species** | Human only | Human + mouse | Human + mouse | Multiple atlases | Human + mouse |
| **Cell type coverage** | 36 common types | Broad ontology | ~400 types | Atlas-dependent | ~500 types |
| **Annotation accuracy** | Not benchmarked at scale | Published benchmarks | Published benchmarks | High (reference-based) | High (reference-based) |

**Where other tools are stronger than MarkerAtlas:**

- **Azimuth** is substantially more accurate for automated annotation of large datasets — it projects your cells onto a curated reference atlas using expression profiles rather than marker overlap, and is the appropriate tool for bulk annotation of standard tissue types (PBMC, brain, lung, kidney, etc.)
- **CellTypist** has broader species coverage (human + mouse), covers ~500 cell types with trained logistic regression models, and has been formally benchmarked
- **sctype** covers ~400 human and mouse cell types from its built-in database, substantially broader than MarkerAtlas's 36
- **ACT** uses a knowledge graph approach that encodes ontological relationships between cell types, which can handle novel or rare cell types better than overlap scoring
- **All of the above** offer or integrate mouse marker support; MarkerAtlas is human-only

**MarkerAtlas is designed for a different use case:** transparent manual annotation with full auditability, where you need to understand and justify each call rather than accept a score, and where you want the annotation workflow — from lookup through batch annotation to figure code — in a single offline tool.

---

## Limitations

**Species:** Human only. Mouse markers differ significantly for some types (Ly6C for monocytes, Siglec-H for pDCs, etc.).

**Coverage:** 36 common cell types. Rare, highly tissue-specific, or recently characterised subtypes are not included (tuft cells, ionocytes, Cajal-Retzius neurons, MAIT cells, etc.).

**Cluster annotation scoring:** Scores marker gene list overlap only — does not account for expression levels, co-expression patterns, or trajectory context. Two clusters with identical top genes score identically regardless of the magnitude of expression differences.

**Evidence tiers:** Based on curated estimates of literature volume, not automated systematic review of publications. Treat as approximate.

**Dynamic and activation-state markers:** Database reflects resting-state consensus annotations. Several markers change meaningfully with activation: IL2RA on activated T cells vs Tregs; GZMB on effector vs exhausted CD8+ T cells; TREM2 on homeostatic vs disease-associated microglia.

**Dropout:** Absence of a marker in sequencing data may reflect technical dropout rather than true biological absence, particularly for lowly expressed genes. Do not use negative marker absence as a hard filter.

**DotPlot code:** Cluster numbering in generated code corresponds to the order of annotation results, not your actual Seurat/Scanpy cluster IDs. Adjust before running.

**PMID hover cards:** Require an internet connection to fetch paper metadata from NCBI. Core annotation functionality works fully offline without them.

---

## Disclaimer

MarkerAtlas is a research reference tool. Marker gene associations are curated from published databases and primary literature but may not reflect all biological contexts, species variants, or disease states. Cell type classifications are based on consensus annotations that may differ from study-specific definitions. Users should verify marker suitability for their specific tissue, species, and experimental conditions before applying annotations in published work. Not intended for clinical, diagnostic, or patient care use.

When using marker information in published work, please cite the primary databases:

> Zhang X et al. CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data. *Nucleic Acids Research*. 2023;51(D1):D870–D876. [PMID 36300619](https://pubmed.ncbi.nlm.nih.gov/36300619/)

> Franzén O, Gan LM, Björkegren JLM. PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. *Database*. 2019;2019:baz046. [PMID 30951143](https://pubmed.ncbi.nlm.nih.gov/30951143/)

---

## License

MIT — free to use, modify, and redistribute with attribution.

---

*Built by [Chakit Arora](https://chakitarora.github.io) · Post-Doctoral Fellow, BIO@SNS, Scuola Normale Superiore, Pisa*

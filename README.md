# MarkerAtlas — Single-Cell Marker Gene Reference

A browser-based reference tool for single-cell RNA sequencing researchers. Five modes: browse all cell types, look up a gene, annotate a cluster, build a minimal marker panel, and inspect annotation ambiguity between cell types. All with verified sources and PubMed links.

**Live tool:** `https://chakitarora.github.io/markeratlas`

Works in any modern browser. Requires a local or hosted server (GitHub Pages, or `python3 -m http.server 8080`) — does not work when opened as a local `file://` URL because `marker_db.json` is loaded via `fetch()`. No account. No data leaves your machine.

---

## The gap this fills

Cell type annotation in scRNA-seq is one of the most consequential and least standardised steps in the workflow. The standard approach — run differential expression, Google the top genes, cross-reference CellMarker or PanglaoDB one gene at a time — is fragmented, slow, and produces undocumented decisions.

Existing tools either require a Python/R environment (CellTypist, SingleR, SCSA), send your data to a remote server (ACT, Azimuth), or show you a score without explaining what matched (ACT, sctype). MarkerAtlas is designed as a transparent offline reference: you see exactly which genes matched, why a cell type was ranked where it was, what you should verify is absent, and which paper supports every association.

It is not a replacement for reference-based automated annotation at scale. It is the tool you reach for when you are looking at a cluster and need to make and justify a biological decision.

---

## Repository contents

```
markeratlas/
├── index.html                  Complete tool — 910 lines, runs on any HTTP server
├── marker_db.json              Full database as JSON (machine-readable, 78 KB)
├── marker_associations.tsv     339 gene–celltype associations as flat TSV
├── celltypes.tsv               36 cell type definitions as flat TSV
└── README.md                   This file
```

---

## How to run locally

```bash
cd markeratlas/
python3 -m http.server 8080
# then open http://localhost:8080 in your browser
```

---

## How to use the tool

The tool has five modes selectable via the toggle buttons below the search box.

---

### Mode 1 — Browse

**What it does:** Displays all 36 cell types as cards. Explore the marker landscape, understand canonical vs supporting vs negative distinctions, and navigate by category or tissue.

**How to use:**
1. Open the tool — Browse is the default
2. Filter by **Category** (Immune, Neural, Stromal, etc.) in the left sidebar
3. Filter by **Tissue** (Blood, Brain, Lung, Tumor, etc.) to show only cell types validated in that context
4. Type in the search box to filter by cell type name or marker gene
5. Click any cell type card to open the full detail panel on the right
6. Click any gene tag on a card to immediately switch to Gene Lookup mode for that gene

**What each card shows:**

| Element | Meaning |
|---|---|
| Blue gene tags | Canonical markers — strong, specific, reproducible evidence |
| Grey gene tags | Supporting markers — useful in combination, not definitive alone |
| Rust/red gene tags | Negative markers — should be absent; important for ruling out contamination |
| Tissue pills | Top 3 validated tissues with confidence % |
| `high` / `medium` badge | Specificity rating — high means reliable across tissues; medium is context-dependent |

**Example — filtering by Neural category + Brain tissue:**

Returns 5 cell types: Neuron, Astrocyte, Oligodendrocyte, Oligodendrocyte precursor cell, Microglia — each with their canonical marker sets and tissue confidence bars.

Clicking **Microglia** opens the detail panel:
```
Minimal panel:  CX3CR1  P2RY12  TMEM119
Canonical:      CX3CR1  P2RY12  TMEM119  HEXB  TREM2
Supporting:     SALL1  FCRLS  ITGAM  CD68
Negative:       GFAP  RBFOX3  OLIG2

Notes: P2RY12 and TMEM119 are highly CNS-specific homeostatic markers,
lost in disease-associated microglia (DAM). TREM2 marks DAM subset.

Tissue confidence: Brain 99%
References: PMID 30726726 ↗  PMID 33558698 ↗
Cell Ontology: CL:0000129 ↗
```

---

### Mode 2 — Gene lookup

**What it does:** Takes a single gene symbol and returns every cell type it marks, with role (canonical or supporting), tissue context, and confidence. Essential for interpreting unexpected DE genes or verifying marker specificity before committing to an annotation.

**How to use:**
1. Click **Gene** in the mode toggle, or click any gene tag in Browse mode
2. Type a gene symbol in the search box — results appear immediately
3. Click any cell type name in the results to open its detail panel

Ten built-in example genes are shown as clickable chips.

**Partial matching:** If the exact symbol is not found, the tool shows genes whose names contain the query string. Typing `SFTPA` returns SFTPB, SFTPC, etc.

**Example — CD3E (specific marker):**

```
CD3E — 3 cell type associations

CD4+ T cell       [canonical]   Blood · Lymph node · Tumor
CD8+ T cell       [canonical]   Blood · Tumor · Lymph node
Regulatory T cell [canonical]   Blood · Tumor · Lymph node
```
CD3E alone identifies T cells but cannot distinguish subtypes. Combine with CD4, CD8A, or FOXP3.

**Example — ACTA2 (ambiguous multi-type marker):**

```
ACTA2 — 4 cell type associations

Smooth muscle cell        [canonical]   Vasculature · Gut · Bladder
Cancer-assoc. fibroblast  [canonical]   Tumor
Fibroblast                [supporting]  Tumor · Skin · Lung
Pericyte                  [supporting]  Brain · Tumor · Kidney
```
ACTA2 (αSMA) marks multiple cell types. Tissue context is essential for correct interpretation. This is exactly the kind of marker that causes annotation errors if treated as specific.

**Example — TREM2 (disease-state marker):**

```
TREM2 — 2 cell type associations

Microglia                [canonical]  Brain
Tumor-assoc. macrophage  [canonical]  Tumor
```
TREM2 marks disease-associated states in both contexts. Its expression indicates an activated/pathological state, not homeostatic identity.

---

### Mode 3 — Cluster annotation

**What it does:** Takes a list of differentially expressed genes from a cluster and scores all 36 cell types by marker overlap. Returns ranked candidates with matched genes, tissue-aware scoring, negative marker suggestions to verify, and an exportable annotation table.

**How to use:**
1. Click **Cluster** in the mode toggle
2. Set the **Tissue context** dropdown to your sample tissue if known (strongly recommended)
3. Paste your top DE genes — one per line, or comma/space separated
4. Click **▶ Annotate cluster**

Six pre-loaded example clusters are available as clickable chips.

**Where to get the input genes:**

```r
# Seurat (R)
markers <- FindMarkers(seurat_obj, ident.1 = 3)
top_genes <- rownames(head(markers[markers$avg_log2FC > 0, ], 20))
cat(top_genes, sep = "\n")
```

```python
# Scanpy (Python)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
top_genes = adata.uns['rank_genes_groups']['names']['3'][:20]
print('\n'.join(top_genes))
```

**Tissue-aware scoring:** When a tissue is selected, each cell type's raw overlap score is multiplied by its tissue confidence for that tissue (0–1). Cell types not documented in that tissue are penalised (×0.1). This means a cluster from a brain sample will not rank Monocytes highly even if CD68 appears — CD68 is a low-weight supporting marker for Microglia, and Monocytes score low in the Brain prior.

Each result shows:
- Rank number (coloured for rank 1)
- Cell type name (clickable → detail panel)
- Tissue prior badge (green = validated in tissue; rust = atypical)
- Match summary (n matched / total markers in database)
- Weighted score
- Matched gene tags
- **Negative marker check** (rust panel): genes that should be absent but were not in your input list — verify these are not expressed in your cluster

**Example — CD8+ T cell cluster, tissue = Blood:**

Input:
```
CD3E, CD3D, CD8A, CD8B, GZMB, GZMK, PRF1, NKG7, GNLY, IFNG
```

Output:
```
1. CD8+ T cell     Blood 99%   10/10 markers · score 18.0  ██████████
   CD3E CD3D CD8A CD8B GZMB GZMK PRF1 NKG7 GNLY IFNG
   → Rule out: verify CD4, CD19, CD14 are absent

2. NK cell         Blood 98%   4/9 markers  · score 6.9    ████
   GZMB PRF1 NKG7 GNLY
   → Rule out: verify CD3E, CD19, CD14 are absent

3. CD4+ T cell     Blood 99%   2/6 markers  · score 4.0    ██
   CD3E CD3D
```

The NK cell at rank 2 is expected — GZMB/PRF1/NKG7/GNLY are shared between cytotoxic T cells and NK cells. CD8A/CD8B are T cell-specific and only appear at rank 1.

**Example — Macrophage cluster, tissue = Tumor:**

Input:
```
C1QA, C1QB, C1QC, CD68, MRC1, APOE, SPP1, TREM2, MARCO, CD163
```

Output:
```
1. Macrophage               Tumor 95%   8/9 · score 12.4   ██████████
   C1QA C1QB C1QC CD68 MRC1 APOE SPP1 TREM2 MARCO
   → Rule out: verify CD3E, CD19, FCGR3B are absent

2. Tumor-assoc. macrophage  Tumor 95%   3/9 · score 5.7    █████
   SPP1 TREM2 CD68

3. Monocyte                 Tumor —     1/8 · score 0.1    █
   CD68  [penalised: not typical in Tumor]
```

SPP1+TREM2 co-expression at rank 2 suggests this may be a disease-associated or tumor-associated macrophage subset. Consider splitting the cluster or checking sample origin.

**Export annotation table:** After annotating, a download button appears. Exports a TSV with rank, cell type, score, matched genes, total markers, tissue prior, and PMIDs — ready to paste into a supplementary table.

---

### Mode 4 — Panel builder

**What it does:** For any cell type, returns the minimal set of genes needed to uniquely identify it, the key negative markers to verify, and a disambiguation guide showing how to distinguish it from the most easily confused neighbours.

Useful for: designing flow cytometry panels, writing methods sections, and deciding which markers to validate by IHC or smFISH.

**How to use:**
1. Click **Panel** in the mode toggle
2. Select a cell type from the grouped dropdown
3. Click **Build panel**
4. Click **↓ Export panel** to download a `.txt` summary

**What the output shows:**

**Minimal panel** — 2–5 genes labelled `unique` (canonical only in this cell type across the full database) or `semi-unique` (canonical in ≤1 other type). These are the genes that carry the most discriminating information. Click any gene to run a gene lookup.

**Negative markers** — genes that should be absent. Shown as rust-coloured tags. Do not include these in a positive staining panel — use them as exclusion criteria.

**Disambiguation table** — the top 3 most easily confused cell types, with:
- Shared markers (the source of potential confusion)
- Positive discriminators (+) — genes in your target type but not the confused type
- Negative discriminators (−) — genes that are negative in your target type but positive in the confused type

**Example — Regulatory T cell:**

```
Minimal panel:
  FOXP3 [unique]   IL2RA [unique]   IKZF2 [unique]   CTLA4 [semi-unique]

Negative markers: CD8A  IL7R

Disambiguation:
  vs CD4+ T cell    shared: CD3E, CD4
                    + FOXP3, IL2RA, IKZF2 (positive discriminators)
                    − IL7R (negative in Treg, positive in CD4+ T cell)

  vs Exhausted T cell   shared: TIGIT, CTLA4
                         + FOXP3, IL2RA (positive discriminators)

  vs CD8+ T cell    shared: CD3E
                    − CD8A (negative in Treg, canonical in CD8+ T cell)
```

**Example — Alveolar type II cell:**

```
Minimal panel:
  SFTPC [unique]   SFTPB [unique]   SFTPD [unique]

Negative markers: AGER  PDPN  TP63

Disambiguation:
  vs Alveolar type I cell   shared: (none in database)
                             + SFTPC, SFTPB, SFTPD
                             − AGER, PDPN (negative in AT2, canonical in AT1)

  vs Club cell   shared: (none in database)
                  + SFTPC, SFTPB
```

All panel genes are clickable and link to the gene lookup view.

---

### Mode 5 — Overlap

**What it does:** Shows the top 20 most ambiguous cell type pairs ranked by number of shared markers. Tells you at a glance which annotation calls are inherently harder to make with marker genes alone, and which genes are responsible for the overlap.

No input required — the table renders immediately when you enter this mode.

**How to use:**
1. Click **Overlap** in the mode toggle
2. Browse the ranked table
3. Click any cell type name to open its detail panel
4. Click any shared gene tag to run a gene lookup

**What the table shows:**

| Column | Meaning |
|---|---|
| Cell type 1 / Cell type 2 | Coloured by category — click to open detail |
| Shared markers | The genes responsible for the overlap — click any to look up |
| N bar | Number of shared markers, shown as a bar chart |

**Example output (top rows):**

```
Macrophage          vs  Tumor-assoc. macrophage   shared: C1QA CD68 MRC1 TREM2 ...  8
CD8+ T cell         vs  NK cell                    shared: GZMB PRF1 NKG7 GNLY      4
Fibroblast          vs  Cancer-assoc. fibroblast   shared: THY1 ACTA2 FAP POSTN     4
Epithelial cell     vs  Cholangiocyte              shared: KRT7 KRT19 EPCAM          3
CD4+ T cell         vs  CD8+ T cell                shared: CD3D CD3E                 2
Regulatory T cell   vs  Exhausted T cell           shared: TIGIT CTLA4               2
Macrophage          vs  Microglia                  shared: TREM2 CD68                 2
Pericyte            vs  Smooth muscle cell         shared: DES ACTA2                  2
```

This tells you directly: macrophage vs TAM is the hardest call in this database (8 shared markers). CD8+ T cell vs NK cell is the classic cytotoxic cell ambiguity. Knowing this going in saves annotation errors.

---

## Detail panel

Clicking any cell type name from any mode opens a sliding detail panel with the complete entry:

| Section | Content |
|---|---|
| **Minimal panel** | 2–5 recommended genes computed from uniqueness across the database |
| **Canonical markers** | Blue clickable tags — each runs a gene lookup |
| **Supporting markers** | Grey clickable tags — each runs a gene lookup |
| **Negative markers** | Rust tags — should be absent |
| **Notes** | Interpretive text: caveats, subtypes, co-expression requirements, disease-state changes |
| **Tissue confidence** | Bar chart: each validated tissue with a 0–100% confidence score |
| **Primary references** | Direct PubMed links (verifiable) |
| **Data sources** | CellMarker 2.0 and/or PanglaoDB |
| **Cell Ontology ID** | CL: accession with link to OLS browser |

---

## Database structure

### Cell types (36 total)

| Category | Cell types |
|---|---|
| Immune (11) | CD4+ T cell, CD8+ T cell, Regulatory T cell, NK cell, B cell, Plasma cell, Monocyte, Macrophage, Dendritic cell, Mast cell, Neutrophil |
| Epithelial (8) | Epithelial cell, Club cell, Alveolar type I, Alveolar type II, Hepatocyte, Cholangiocyte, Enterocyte, Goblet cell |
| Stromal (5) | Fibroblast, Endothelial cell, Pericyte, Smooth muscle cell, Adipocyte |
| Neural (5) | Neuron, Astrocyte, Oligodendrocyte, OPC, Microglia |
| Stem/Progenitor (2) | Intestinal stem cell, Hematopoietic stem cell |
| Tumor microenvironment (3) | Cancer-associated fibroblast, Tumor-associated macrophage, Exhausted T cell |
| Hematopoietic (2) | Erythrocyte, Megakaryocyte |

### Marker confidence tiers

**Canonical** — validated in multiple independent studies, listed in both CellMarker and PanglaoDB, with primary PMIDs. These define the cell type.

**Supporting** — useful in combination but insufficient alone. May mark multiple cell types or have tissue-specific expression.

**Negative** — absence is diagnostically important. Dropout in sequencing data may mimic true biological absence — use negative markers as a verification check, not a hard filter.

### Database fields (marker_db.json)

Each cell type entry contains:
```json
{
  "category": "Immune",
  "canonical_markers": ["CD3E", "CD3D", "CD8A", ...],
  "supporting_markers": ["PRF1", "NKG7", ...],
  "negative_markers": ["CD4", "CD19", ...],
  "tissue_context": {"Blood": 0.99, "Tumor": 0.95, ...},
  "specificity": "high",
  "notes": "...",
  "pmids": ["30726726", ...],
  "sources": ["CellMarker", "PanglaoDB"],
  "cl_id": "CL:0000625",
  "minimal_panel": ["CD8A", "CD8B", "GZMK"],
  "unique_markers": ["CD8A", "CD8B"]
}
```

`minimal_panel` — 2–5 genes combining unique and semi-unique canonical markers.
`unique_markers` — canonical markers that appear in no other cell type's canonical list.
`cl_id` — Cell Ontology accession, linkable to OLS4 browser.

---

## Using the data files in R / Python

```r
# R
markers <- read.delim("marker_associations.tsv", stringsAsFactors=FALSE)

# All canonical markers for immune types
subset(markers, category=="Immune" & role=="canonical")

# All cell types for a given gene
subset(markers, gene=="TREM2")

# Genes unique to one cell type
subset(markers, specificity=="high" & role=="canonical")[
  !duplicated(subset(markers, specificity=="high" & role=="canonical")$gene), ]
```

```python
# Python
import pandas as pd, json

# TSV approach
m = pd.read_csv("marker_associations.tsv", sep="\t")
m[m.gene == "FOXP3"]  # all associations for FOXP3

# JSON approach — access minimal panels
with open("marker_db.json") as f:
    db = json.load(f)
panels = {ct: data["minimal_panel"] for ct, data in db["celltypes"].items()}
```

---

## Data sources

| Database | Reference | URL |
|---|---|---|
| CellMarker 2.0 | Zhang X et al., *Nucleic Acids Res.* 2023. [PMID 36300619](https://pubmed.ncbi.nlm.nih.gov/36300619/) | [biocc.hrbmu.edu.cn/CellMarker](https://biocc.hrbmu.edu.cn/CellMarker/) |
| PanglaoDB | Franzén O et al., *Database* 2019. [PMID 30951143](https://pubmed.ncbi.nlm.nih.gov/30951143/) | [panglaodb.se](https://panglaodb.se) |

---

## Comparison with related tools

| | MarkerAtlas | ACT | sctype | Azimuth | CellTypist |
|---|---|---|---|---|---|
| Input | Marker gene list | Marker gene list | Marker gene list | Expression matrix | Expression matrix |
| Infrastructure | Browser, offline | Remote server | Browser / R | Remote server | Python package |
| Transparency | Full — shows matched genes | Score only | Partial | Score only | Score only |
| Negative markers | Yes | No | No | No | No |
| Tissue-aware scoring | Yes | No | No | Yes | No |
| Minimal panel builder | Yes | No | No | No | No |
| Confusion matrix | Yes | No | No | No | No |
| Verifiable PMIDs | Yes | No | No | No | No |
| Offline operation | Yes | No | No | No | No |
| Species | Human | Human + mouse | Human + mouse | Multiple atlases | Human + mouse |

MarkerAtlas is a **transparent reference tool for manual annotation**, not a replacement for reference-based automated annotation at scale.

---

## Limitations

**Species** — human (*Homo sapiens*) only. Mouse markers differ for some types (Ly6C for monocytes, Tmem119 capitalisation, etc.).

**Coverage** — 36 common cell types. Rare subtypes (tuft cells, ionocytes, Cajal-Retzius cells) not included.

**Cluster annotation scoring** — scores marker overlap only. Does not account for expression levels, co-expression, or trajectory. Use as a starting hypothesis.

**Dynamic markers** — some markers change with activation state (IL2RA on activated vs Treg; GZMB on effector vs exhausted CD8+). Database reflects resting-state consensus.

**Dropout** — absence of a marker in sequencing data may be technical. Do not use negative markers as hard exclusion filters.

---

## Disclaimer

MarkerAtlas is a research reference tool. Marker gene associations are curated from published databases and primary literature but may not reflect all biological contexts, species variants, or disease states. Cell type classifications are based on consensus annotations that may differ from study-specific definitions. Users should verify marker suitability for their specific tissue, species, and experimental conditions before applying annotations in published work. Not intended for clinical, diagnostic, or patient care use.

When using marker information in published work, cite the primary databases (CellMarker 2.0, PanglaoDB) and the listed PMIDs.

---

## Deploy

```bash
# GitHub Pages
# Fork or upload this repo, go to Settings → Pages → Branch: main → Save
# Live at https://yourusername.github.io/markeratlas
```

Runs entirely from `index.html` + `marker_db.json`. No build step. No server-side logic.

---

## License

MIT — free to use, modify, and redistribute with attribution.

---

*Built by [Chakit Arora](https://chakitarora.github.io) · Post-Doctoral Fellow, BIO@SNS, Scuola Normale Superiore, Pisa*

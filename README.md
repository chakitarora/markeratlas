# MarkerAtlas — Single-Cell Marker Gene Reference

A browser-based reference tool for single-cell RNA sequencing researchers. Search a gene, browse cell types, or paste a cluster's top DE genes to get ranked candidate annotations — all with verified sources and PubMed links.

**Live tool:** `https://chakitarora.github.io/markeratlas`

Works in any modern browser. No internet required after the page loads. No account. No data leaves your machine.

---

## The gap this fills

Cell type annotation is one of the most consequential and least standardised steps in scRNA-seq analysis. The standard workflow — run differential expression, Google the top genes, cross-reference CellMarker or PanglaoDB manually — is fragmented, slow, and produces undocumented annotation decisions.

Existing databases (CellMarker 2.0, PanglaoDB) are excellent resources but are designed for one-gene-at-a-time lookup. Neither can answer: *"Given these 15 DE genes from my cluster, what cell type am I looking at, and how confident should I be?"* They also do not surface conflicts — genes that mark multiple cell types where tissue context determines the correct interpretation.

MarkerAtlas aggregates both databases into a single interface with three complementary search modes, confidence scoring, and explicit negative markers.

---

## Repository contents

```
markeratlas/
├── index.html                  Complete tool — self-contained, works offline
├── marker_db.json              Full database as JSON (machine-readable)
├── marker_associations.tsv     339 gene–celltype associations as flat TSV
├── celltypes.tsv               36 cell type definitions as flat TSV
└── README.md                   This file
```

---

## How to use the tool

The tool has three modes, selectable via the toggle buttons below the search box.

---

### Mode 1 — Browse

**What it does:** Displays all 36 cell types as cards. Use it to explore what markers define a cell type, understand canonical vs supporting vs negative marker distinctions, and navigate to cell types by category or tissue.

**How to use:**
1. Open the tool — Browse mode is the default
2. Use the **Category** filter in the left sidebar to narrow to a group (Immune, Neural, Stromal, etc.)
3. Use the **Tissue** filter to show only cell types validated in a specific tissue (Blood, Brain, Lung, Tumor, etc.)
4. Type in the search box to filter by cell type name or marker gene
5. Click any cell type card to open the full detail panel on the right

**What each card shows:**

| Element | Meaning |
|---|---|
| Blue gene tags | Canonical markers — strong, specific evidence. Use these to define the cell type. |
| Grey gene tags | Supporting markers — useful in combination but not definitive alone |
| Rust/red gene tags | Negative markers — should be absent. Important for ruling out contamination or doublets. |
| Tissue pills | Top 3 tissues where validated, with confidence % |
| Specificity badge | `high` = reliable across tissues; `medium` = context-dependent |

Clicking any **gene tag** immediately switches to Gene Lookup mode for that gene.

**Example — searching "T cell" in Browse mode:**

Typing `T cell` filters to: CD4+ T cell, CD8+ T cell, Regulatory T cell, Exhausted T cell.

Clicking **CD8+ T cell** opens the detail panel:
```
Canonical:   CD3E  CD3D  CD8A  CD8B  GZMK  GZMB
Supporting:  PRF1  NKG7  GNLY  IFNG
Negative:    CD4   CD19  CD14

Notes: Cytotoxic T cells. GZMB/PRF1 mark effector function. Exhausted
CD8+ T cells additionally express PDCD1, HAVCR2, TIGIT.

Tissue confidence:
  Blood       99%  █████████░
  Tumor       95%  █████████░
  Lymph node  97%  █████████░

References: PMID 30726726 ↗  PMID 33558698 ↗  PMID 32103181 ↗
```

---

### Mode 2 — Gene lookup

**What it does:** Takes a single gene symbol and returns every cell type it has been reported to mark, with role (canonical or supporting), tissue context, and confidence. Essential for interpreting unexpected DE genes or checking marker specificity.

**How to use:**
1. Click **Gene lookup** in the mode toggle
2. Type a gene symbol in the search box
3. Results appear immediately
4. Click any cell type name to open its full detail panel

Built-in example genes are shown as clickable chips — CD3E, MBP, SFTPC, FOXP3, SPP1, PDGFRB, LGR5, TREM2, ACTA2, NCAM1.

**Example 1 — specific marker (CD3E):**

Input: `CD3E`

```
CD3E — 3 cell type associations

CD4+ T cell       [canonical]   Blood · Lymph node · Tumor
CD8+ T cell       [canonical]   Blood · Tumor · Lymph node
Regulatory T cell [canonical]   Blood · Tumor · Lymph node
```

Interpretation: CD3E is a pan-T cell marker. Alone it cannot distinguish CD4+, CD8+, or Treg. Always combine with CD4, CD8A, or FOXP3 respectively.

**Example 2 — multi-type marker (ACTA2):**

Input: `ACTA2`

```
ACTA2 — 4 cell type associations

Smooth muscle cell         [canonical]   Vasculature · Gut · Bladder
Cancer-assoc. fibroblast   [canonical]   Tumor
Fibroblast                 [supporting]  Tumor · Skin · Lung
Pericyte                   [supporting]  Brain · Tumor · Kidney
```

Interpretation: ACTA2 (αSMA) marks multiple cell types. In vasculature it marks smooth muscle; in tumor stroma it marks CAFs; in brain it marks pericytes. Tissue context is essential.

**Example 3 — disease-state marker (TREM2):**

Input: `TREM2`

```
TREM2 — 2 cell type associations

Microglia                  [canonical]  Brain
Tumor-assoc. macrophage    [canonical]  Tumor
```

Interpretation: TREM2 marks disease-associated microglia (DAM) in brain and immunosuppressive TAMs in tumors. Its expression indicates an activated/pathological state in both contexts.

**Partial matching:** If the exact symbol is not found, the tool shows partial matches. Typing `CD8` returns CD8A, CD8B.

---

### Mode 3 — Cluster annotation

**What it does:** Takes a list of differentially expressed genes from a single cluster and scores all 36 cell types by marker overlap. Returns ranked candidate annotations showing which of your genes matched each cell type, weighted by marker role (canonical = 2pts, supporting = 1pt).

**How to use:**
1. Click **Cluster** in the mode toggle
2. Paste your top DE genes — one per line, or comma/space separated
3. Click **▶ Annotate cluster**
4. Ranked results appear with scores, match summaries, and matched gene tags

Six pre-loaded example clusters are available as clickable buttons: CD8+ T cell, Macrophage, Astrocyte, AT2 lung cell, Fibroblast, and an Ambiguous cluster.

**Where to get the input genes:**

From Seurat (R):
```r
markers <- FindMarkers(seurat_obj, ident.1 = 3)
top_genes <- rownames(head(markers[markers$avg_log2FC > 0, ], 20))
cat(top_genes, sep = "\n")
```

From Scanpy (Python):
```python
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
top_genes = adata.uns['rank_genes_groups']['names']['3'][:20]
print('\n'.join(top_genes))
```

**Example 1 — clear annotation (CD8+ T cell cluster):**

Input:
```
CD3E
CD3D
CD8A
CD8B
GZMB
GZMK
PRF1
NKG7
GNLY
IFNG
```

Output:
```
1. CD8+ T cell         10/10 markers · score 18  ██████████
   CD3E CD3D CD8A CD8B GZMB GZMK PRF1 NKG7 GNLY IFNG

2. NK cell             4/9 markers  · score 7    ████
   GZMB PRF1 NKG7 GNLY

3. CD4+ T cell         2/6 markers  · score 4    ██
   CD3E CD3D

4. Regulatory T cell   1/6 markers  · score 1    █
   CD3E
```

Interpretation: High confidence CD8+ T cell. NK cell at rank 2 is expected — GZMB/PRF1/NKG7/GNLY are shared between cytotoxic T cells and NK cells, but CD8A/CD8B are T cell-specific and only appear at rank 1.

**Example 2 — subtype context (Macrophage cluster):**

Input:
```
C1QA
C1QB
C1QC
CD68
MRC1
APOE
SPP1
TREM2
MARCO
CD163
```

Output:
```
1. Macrophage              8/9 markers · score 13  ██████████
   C1QA C1QB C1QC CD68 MRC1 APOE SPP1 TREM2 MARCO

2. Tumor-assoc. macrophage 3/9 markers · score 6   █████
   SPP1 TREM2 CD68

3. Monocyte                1/8 markers · score 1   █
   CD68
```

Interpretation: Macrophage with high confidence. The SPP1 + TREM2 co-expression at rank 2 suggests tumor-associated or disease-associated macrophage subset — consider the sample origin. If this is from tumor tissue, this cluster may warrant splitting.

**Example 3 — ambiguous/doublet cluster:**

Input:
```
CD3E
COL1A1
GFAP
EPCAM
CD68
PECAM1
HBB
SFTPC
```

Output:
```
1. Fibroblast      1/9 markers · score 1
2. Epithelial cell 1/7 markers · score 1
3. Macrophage      1/9 markers · score 1
4. Hepatocyte      1/8 markers · score 1
   (all scores equal and low)
```

Interpretation: No dominant cell type — genes span immune, stromal, epithelial, and neural lineages simultaneously. This is biologically implausible for a single population and strongly suggests doublets or a contaminated cluster. Run DoubletFinder (R) or Scrublet (Python) and re-cluster.

---

## Detail panel

Clicking any cell type name from any mode opens a sliding detail panel:

| Section | Content |
|---|---|
| Canonical markers | Blue clickable tags — each runs a gene lookup |
| Supporting markers | Grey clickable tags — each runs a gene lookup |
| Negative markers | Rust tags — genes that should be absent |
| Notes | Caveats, subtypes, co-expression requirements, disease-state changes |
| Tissue confidence | Bar chart with 0–100% confidence per tissue |
| Primary references | Direct PubMed links (verifiable) |
| Data sources | CellMarker 2.0 and/or PanglaoDB |

---

## Database coverage

**36 cell types · 289 unique marker genes · 339 gene–celltype associations**

All from human (*Homo sapiens*). Sources: CellMarker 2.0 (Zhang et al. 2023) and PanglaoDB (Franzén et al. 2019), with primary literature PMIDs per entry.

Full cell type list by category:

| Category | Cell types |
|---|---|
| Immune | CD4+ T cell, CD8+ T cell, Regulatory T cell, NK cell, B cell, Plasma cell, Monocyte, Macrophage, Dendritic cell, Mast cell, Neutrophil |
| Epithelial | Epithelial cell, Club cell, Alveolar type I, Alveolar type II, Hepatocyte, Cholangiocyte, Enterocyte, Goblet cell |
| Stromal | Fibroblast, Endothelial cell, Pericyte, Smooth muscle cell, Adipocyte |
| Neural | Neuron, Astrocyte, Oligodendrocyte, Oligodendrocyte precursor cell, Microglia |
| Stem/Progenitor | Intestinal stem cell, Hematopoietic stem cell |
| Tumor microenvironment | Cancer-associated fibroblast, Tumor-associated macrophage, Exhausted T cell |
| Hematopoietic | Erythrocyte, Megakaryocyte |

---

## Using the data files in R / Python

```r
# R — load marker associations
markers <- read.delim("marker_associations.tsv", stringsAsFactors=FALSE)

# All canonical markers for immune cell types
subset(markers, category=="Immune" & role=="canonical")

# Which cell types does PDGFRB mark?
subset(markers, gene=="PDGFRB")

# All high-specificity markers
subset(markers, specificity=="high" & role=="canonical")
```

```python
# Python — load with pandas
import pandas as pd
markers = pd.read_csv("marker_associations.tsv", sep="\t")

# Cell types for a list of genes from a cluster
cluster_genes = ["CD3E","CD8A","GZMB","PRF1","NKG7"]
markers[markers.gene.isin(cluster_genes)].groupby("celltype").size().sort_values(ascending=False)
```

---

## Limitations

**Species:** Human only. Mouse markers differ for some types.

**Coverage:** 36 common cell types. Rare/tissue-specific subtypes (tuft cells, ionocytes, Cajal-Retzius cells) not included.

**Cluster annotation:** Scores marker overlap only — does not account for expression levels, co-expression, or trajectory. A starting hypothesis, not a definitive annotation.

**Dynamic markers:** Some markers change with activation (IL2RA on activated vs Treg; GZMB on effector vs exhausted CD8+). Database reflects resting-state consensus.

**Dropout:** Absence of a marker in sequencing data may reflect technical dropout rather than true biological absence, especially for lowly expressed genes.

---

## Disclaimer

MarkerAtlas is a research reference tool. Marker associations are curated from published databases and primary literature but may not reflect all biological contexts, species variants, or disease states. Cell type classifications are based on consensus annotations and may differ from study-specific definitions. Users should verify marker suitability for their specific tissue, species, and experimental conditions before applying annotations in published work. This tool is not intended for clinical, diagnostic, or patient care use.

---

## Data sources

| Database | Reference | URL |
|---|---|---|
| CellMarker 2.0 | Zhang X et al., *Nucleic Acids Res.* 2023. [PMID 36300619](https://pubmed.ncbi.nlm.nih.gov/36300619/) | [biocc.hrbmu.edu.cn/CellMarker](https://biocc.hrbmu.edu.cn/CellMarker/) |
| PanglaoDB | Franzén O et al., *Database* 2019. [PMID 30951143](https://pubmed.ncbi.nlm.nih.gov/30951143/) | [panglaodb.se](https://panglaodb.se) |

---

## License

MIT — free to use, modify, and redistribute with attribution.

---

*Built by [Chakit Arora](https://chakitarora.github.io) · Post-Doctoral Fellow, BIO@SNS, Scuola Normale Superiore, Pisa*

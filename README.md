# MarkerAtlas — Single-Cell Marker Gene Reference

A browser-based reference tool for single-cell RNA-seq cell type annotation. Nine modes covering the full annotation workflow: browsing the marker database, gene-level lookups, cluster annotation, batch annotation, co-expression queries, minimal panel design, ambiguity analysis, dot plot code generation, and cell type name standardisation against the Cell Ontology.

**Live tool:** `https://chakitarora.github.io/markeratlas`

Requires `marker_db.json` in the same directory. Works on GitHub Pages or locally with `python3 -m http.server 8080`. No account. No install. No data leaves your browser.

---

## What problem this solves

Annotating clusters in scRNA-seq analysis involves three recurring bottlenecks:

1. **Which cell type do these marker genes point to?** — Cluster annotation mode scores your DE gene list against the database.
2. **Is this gene a good marker, and for what?** — Gene mode shows every cell type a gene marks, its tier, and whether it is specific or shared.
3. **My collaborators and I used different names for the same cell type** — Standardise mode maps any informal label (Mφ, TAM, Treg, CD56+ lymphocyte) to the Cell Ontology standard name and CL ID, detects duplicates, and generates renaming code.

---

## Database

**36 cell types** across 7 categories: Immune (T, B, NK, myeloid), Stromal, Epithelial, Neural, Endothelial, Progenitor, Other.

**289 genes** in the gene index with evidence tiers.

**Fields per cell type:**
- `canonical_markers` — gold-standard identifying genes
- `supporting_markers` — supporting but not exclusive
- `negative_markers` — genes that should be absent
- `tissue_context` — per-tissue confidence scores
- `specificity` — how unique the marker set is
- `cl_id` — Cell Ontology identifier
- `minimal_panel` — smallest gene set that uniquely identifies this type
- `unique_markers` — genes exclusive to this cell type
- `pmids` — PubMed references supporting the annotation

**Gene evidence tiers:**
- `gold` — canonical, widely validated, high specificity
- `silver` — well-supported, moderate specificity
- `bronze` — supporting evidence, lower specificity or tissue-restricted

---

## The nine modes

---

### 1 — Browse

**When to use:** exploring the database, checking what markers are associated with a cell type, or finding a starting point for a new tissue.

Displays all 36 cell types as cards. Filter by **category** (Immune, Stromal, Epithelial, etc.) and **tissue** using the sidebar controls.

Each card shows: cell type name, category badge, canonical marker genes (colour-coded by tier), tissue contexts, and specificity rating.

**Click any gene tag** to run a gene lookup for that marker.
**Click a card** to open the full detail panel showing:
- All marker tiers with evidence badges
- Tissue confidence scores with confidence bars
- Negative markers
- Minimal panel genes
- PubMed references (hover for title)
- CL ID with link to OLS4

---

### 2 — Gene

**When to use:** you have a gene and want to know which cell types it marks and how specifically.

Type any gene symbol. **Aliases are resolved** — the following informal and CD nomenclature names work:

| Alias | Resolves to |
|---|---|
| CD45 | PTPRC |
| PD-1 | PDCD1 |
| PD-L1 | CD274 |
| CD56 | NCAM1 |
| SMA | ACTA2 |
| NeuN | RBFOX3 |
| Iba1 | AIF1 |
| GFAP | GFAP |
| Ki67 | MKI67 |
| EpCAM | EPCAM |

Results show every cell type in the database that lists the gene, with:
- Role: **canonical**, **supporting**, or **negative**
- Evidence tier badge (gold / silver / bronze)
- Tissue context
- Cell type specificity

Results sorted by tier then specificity. Click any cell type name to open its detail panel.

**Partial match:** typing 3+ characters shows partial matches before exact results.

---

### 3 — Cluster annotation

**When to use:** annotating one cluster at a time from your DE gene list.

Paste your top differentially expressed genes — one per line or comma-separated. Select a tissue context from the dropdown to apply tissue-specific priors (improves accuracy for tissue-restricted types).

**How scoring works:**

Each cell type receives a score based on:
- Number of canonical markers matched (weighted higher)
- Number of supporting markers matched
- Tissue prior multiplier (if tissue selected and cell type is known for that tissue)
- Penalty for negative markers present in your gene list

Results ranked by score. For each candidate:
- Green ticks show matched canonical markers
- Grey ticks show matched supporting markers
- Red crosses show negative markers found (if any) — these are warnings that the annotation may be wrong

**Confidence thresholds:**
- **High** (≥70% of canonical markers matched, no negative markers)
- **Medium** (≥40% matched, or one negative marker)
- **Low** (≤40% matched, or multiple negative markers)

A threshold slider adjusts the minimum score shown. Session history keeps the last 5 annotations for comparison.

**Export TSV** — one row per result, includes matched genes and confidence.

**Example input (PBMC CD8 T cell cluster):**
```
CD8A
CD8B
GZMB
PRF1
NKG7
IFNG
CD3D
CD3E
```

---

### 4 — Batch cluster annotation

**When to use:** annotating all clusters at once from a full DE analysis output.

**Input format:** one cluster per line, tab-separated:
```
Cluster0    CD8A,CD8B,GZMB,PRF1,NKG7
Cluster1    CD4,IL7R,CCR7,TCF7,LEF1
Cluster2    CD14,LYZ,CST3,FCGR3A,MS4A7
```

Or paste gene-list only rows (clusters numbered automatically):
```
CD8A,CD8B,GZMB,PRF1
CD4,IL7R,CCR7,TCF7
CD14,LYZ,CST3
```

Results table shows: cluster ID, top annotation, confidence rating (High/Medium/Low), matched markers, and top alternative annotation.

**Export TSV** — ready to use as a supplementary table or to import as metadata in Seurat/Scanpy.

---

### 5 — Co-expression query

**When to use:** you want to find cell types where a specific combination of genes is co-expressed — a question cluster annotation alone cannot reliably answer.

Enter 2–6 gene symbols separated by commas or newlines.

**Example:** `CD3E, CD8A, GZMB, PRF1`

Results:
- **Full matches** — cell types where ALL entered genes are listed (canonical or supporting)
- **Partial matches** — cell types where most but not all genes are listed, with missing genes struck through

Each result shows the role (C = canonical, S = supporting) of each gene in that cell type.

**Why this is useful:** two cell types may both annotate as "CD8 T cell" from cluster annotation, but one co-expresses GZMB+PRF1 (cytotoxic effector) and one does not (naive/memory). Co-expression query distinguishes them before you run a sub-clustering step.

---

### 6 — Minimal panel

**When to use:** designing a targeted gene panel, a flow cytometry antibody panel, or writing the methods section of a paper describing how you identified a cell type.

Select a cell type from the dropdown.

**Output:**
- **Minimal panel genes** — the smallest set that uniquely or semi-uniquely identifies this cell type
  - `unique` badge — this gene is listed only for this cell type in the database
  - `semi-unique` badge — listed for at most one other cell type
- **Key negative markers** — genes that should be absent, with the cell types they rule out
- **Disambiguation table** — for the 3 most easily confused neighbour cell types: which genes distinguish this type from each neighbour

**Example (CD8+ T cell):**
```
Minimal panel: CD8A (unique), CD8B (unique), GZMB (semi-unique)
Key negatives: CD4 (rules out CD4+ T), CD19 (rules out B cell)
Disambiguation from NK cell: CD3D present, NKP46 absent
```

**Export panel genes** as a plain text list.

---

### 7 — Overlap (ambiguity analysis)

**When to use:** before starting annotation, to understand which cell types your marker database will struggle to distinguish.

Displays the **20 most ambiguous cell type pairs** ranked by number of shared canonical and supporting markers. For each pair:
- Number of shared genes
- List of shared gene symbols (click any to run a gene lookup)
- Number of unique genes per type (click cell type names to open detail)

**Why this matters:** if Macrophage and Dendritic cell share 8 markers, you should expect annotation ambiguity between them and may need additional discriminating genes in your panel or a sub-clustering step. Knowing this before annotating prevents post-hoc rationalisation.

---

### 8 — DotPlot code

**When to use:** after completing cluster annotation, to generate publication-ready dot plot code.

This mode reads the session state from the most recent **Cluster annotation** or **Batch annotation** run.

Select:
- **Language:** Seurat (R) or Scanpy (Python)
- **Markers per cell type:** 3, 5, or 8

**Output:** complete, runnable code including:
- Cluster renaming using your annotations
- `DotPlot()` call (Seurat) or `sc.pl.dotplot()` call (Scanpy)
- Gene list constructed from top-tier markers for each annotated cell type

Paste directly into your analysis notebook. The gene list is deduplicated and ordered by cell type.

---

### 9 — Standardise

**When to use:** you have cell type names from your own analysis, a collaborator, or a published paper and need to map them to consistent Cell Ontology (CL) standard names before combining datasets, depositing to CELLxGENE, or writing a manuscript.

Paste cell type names one per line — any format works:

```
Mφ
tissue macrophage
NK cell
CD56+ lymphocyte
TAM
CD8 T cell
Treg
activated B
```

**Load example sets:** four preset example datasets are available (PBMC, Tumor microenvironment, Brain atlas, Lung atlas) — click any to pre-fill the input and run.

**Output — three components:**

**1. Standardisation table**

| Input | Standard CL name | CL ID | Confidence | Parent |
|---|---|---|---|---|
| Mφ | Macrophage | CL:0000235 | high | Mononuclear phagocyte |
| tissue macrophage | Macrophage | CL:0000235 | high | Mononuclear phagocyte |
| NK cell | Natural killer cell | CL:0000623 | high | Lymphocyte |
| CD56+ lymphocyte | Natural killer cell | CL:0000623 | medium | Lymphocyte |
| TAM | Tumor-associated macrophage | CL:0001056 | high | Macrophage |
| CD8 T cell | CD8-positive T cell | CL:0000625 | high | T cell |
| Treg | Regulatory T cell | CL:0000815 | high | CD4-positive T cell |
| activated B | — | — | low | ambiguous — activation state |

CL IDs link to the OLS4 Cell Ontology browser.

**Confidence levels:**
- `high` — exact synonym match from curated table (376 synonyms across 36 cell types)
- `medium` — partial string match with one plausible candidate
- `low` — ambiguous or multiple candidates
- `none` — no match found

**2. Hierarchy tree**

Indented tree grouping your cell types by parent lineage:

```
Lymphocyte
  └ Natural killer cell  CL:0000623  ← NK cell, CD56+ lymphocyte
  └ CD4-positive T cell  CL:0000492
    └ Regulatory T cell  CL:0000815  ← Treg
  └ CD8-positive T cell  CL:0000625  ← CD8 T cell

Mononuclear phagocyte
  └ Macrophage           CL:0000235  ← Mφ, tissue macrophage  [duplicate]
    └ Tumor-assoc. mac.  CL:0001056  ← TAM
```

**3. Consistency flags**

- ⚠ "Mφ" and "tissue macrophage" both resolve to Macrophage (CL:0000235) — **duplicate annotation**
- ⚠ "NK cell" and "CD56+ lymphocyte" both resolve to Natural killer cell — **duplicate annotation**
- ⚠ "activated B" — activation state, not a CL lineage term — consider B cell (CL:0000236)
- ℹ 1 input had no CL match

**Two exports:**
- **Export TSV** — full standardisation table with CL IDs, parent terms, and confidence
- **Export R code** or **Export Python code** — ready-to-run renaming script:

```r
# Seurat
cell_type_map <- c(
  "Mφ"                = "Macrophage",
  "tissue macrophage" = "Macrophage",
  "NK cell"           = "Natural killer cell",
  "CD56+ lymphocyte"  = "Natural killer cell",
  "TAM"               = "Tumor-associated macrophage",
  "CD8 T cell"        = "CD8-positive T cell",
  "Treg"              = "Regulatory T cell"
)
seurat_obj$cell_type_standard <- dplyr::recode(seurat_obj$cell_type, !!!cell_type_map)
seurat_obj$cell_type_cl_id    <- dplyr::recode(seurat_obj$cell_type_standard, !!!cl_id_map)
```

```python
# Scanpy
cell_type_map = {
    "Mφ":                "Macrophage",
    "tissue macrophage": "Macrophage",
    "NK cell":           "Natural killer cell",
    "CD56+ lymphocyte":  "Natural killer cell",
    "TAM":               "Tumor-associated macrophage",
    "CD8 T cell":        "CD8-positive T cell",
    "Treg":              "Regulatory T cell",
}
adata.obs["cell_type_standard"] = adata.obs["cell_type"].map(cell_type_map).fillna(adata.obs["cell_type"])
adata.obs["cell_type_cl_id"]    = adata.obs["cell_type_standard"].map(cl_id_map)
```

**Covered cell types (36):** CD4+ T cell, CD8+ T cell, Regulatory T cell, Exhausted T cell, NK cell, B cell, Plasma cell, Monocyte, Macrophage, Tumor-associated macrophage, Dendritic cell, Mast cell, Neutrophil, Fibroblast, Cancer-associated fibroblast, Endothelial cell, Pericyte, Smooth muscle cell, Epithelial cell, Alveolar type I/II cell, Hepatocyte, Cholangiocyte, Enterocyte, Goblet cell, Club cell, Neuron, Astrocyte, Oligodendrocyte, OPC, Microglia, Intestinal stem cell, Hematopoietic stem cell, Erythrocyte, Megakaryocyte, Adipocyte.

---

## Workflow integration

**Typical scRNA-seq annotation workflow with MarkerAtlas:**

```
1. Cluster → FindMarkers (Seurat) or rank_genes_groups (Scanpy)
2. Paste top markers per cluster into Batch annotation
3. Review results — check Overlap mode for ambiguous pairs
4. Use Gene mode for any uncertain markers
5. Use Co-expression query to distinguish subtypes
6. Use Panel mode when writing methods ("CD8+ T cells were identified by...")
7. Use DotPlot mode to generate publication figure code
8. Use Standardise mode to harmonise names before:
   - combining with a collaborator's dataset
   - depositing to CELLxGENE (requires CL IDs)
   - submitting manuscript (consistent nomenclature)
```

---

## Data files

| File | Contents |
|---|---|
| `markeratlas.html` | Complete tool — all modes, all JS, embedded CL database |
| `marker_db.json` | Marker database — 36 cell types, 289 genes, all fields |
| `marker_associations.tsv` | Flat TSV of all gene-cell type associations with tiers |
| `celltypes.tsv` | Flat TSV of cell type metadata including CL IDs |

`markeratlas.html` must be served alongside `marker_db.json` — it fetches the database via HTTP on load. Directly opening the HTML file (`file://`) will fail due to browser CORS restrictions on local file fetches. Use GitHub Pages or `python3 -m http.server 8080`.

---

## Limitations

**36 cell types** — covers major types from blood, brain, gut, lung, liver, and general stromal lineages. Does not cover rare subtypes (e.g. ILC subsets, hair follicle cell types, retinal cell types) or highly tissue-specific types from specialised organs.

**289 marker genes** — well-characterised markers only. Novel markers from recent preprints or highly tissue-specific expression patterns may not be included.

**Standardise mode — 376 synonyms** — covers common informal names and abbreviations in the scRNA-seq literature. Highly context-specific informal names from individual labs may not match. The `none` confidence result is a signal to check the OLS4 browser manually.

**Cluster annotation scoring** — based on marker overlap, not expression level, trajectory, or spatial context. Two clusters with identical top markers will score identically regardless of expression magnitude. Use sub-clustering and additional discriminating markers for fine-grained subtypes.

**Tissue priors** — based on known tissue distributions in the literature. Novel tissue-specific cell types or disease-context populations (e.g. stress-response macrophage states) may score poorly.

---

## Citing MarkerAtlas

When using MarkerAtlas annotations in published work, cite the primary literature from which markers are drawn (PubMed references are shown in each cell type detail panel). For the Cell Ontology standard names:

> Diehl AD et al. The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability. *Journal of Biomedical Semantics.* 2016;7:44.

For the CELLxGENE schema and CL term requirements:

> Megill C et al. CELLxGENE: a performant, scalable exploration of large-scale single-cell data. *bioRxiv.* 2021.

---

## Design

MarkerAtlas uses the Kandinsky Composition VIII visual language: warm canvas (`#f4f0e8`), geometric elements, IBM Plex Mono for identifiers, IBM Plex Sans for body text. Colour encodes evidence tier (gold / silver / bronze) and confidence level throughout.

---

## License

MIT

---

*Built by [Chakit Arora](https://chakitarora.github.io) · Post-Doctoral Fellow, BIO@SNS, Scuola Normale Superiore, Pisa*

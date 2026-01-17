# Transcriptomic Comparison of Chronic Hypersensitivity Pneumonitis and Idiopathic Pulmonary Fibrosis  
### A Reproduction Study of GSE150910 Using DESeq2

## Overview
This repository contains a **reproducible RNA-sequencing differential expression analysis** comparing **Chronic Hypersensitivity Pneumonitis (CHP)**, **Idiopathic Pulmonary Fibrosis (IPF)**, and **healthy control lung tissue** using the publicly available dataset **GSE150910**.

The goal of this work is to **reproduce and validate the key molecular findings** of the original study by identifying:
- Gene expression changes associated with fibrotic lung disease
- Molecular signatures shared between CHP and IPF
- Genes that distinguish CHP from IPF

This project follows a **transparent, step-by-step computational pipeline** using **DESeq2**, emphasizing quality control, interpretability, and reproducibility.

## Scientific Motivation
CHP and IPF are both fibrotic interstitial lung diseases with overlapping clinical and histopathological features, yet they differ in etiology, prognosis, and treatment response. Understanding whether these diseases share a **common fibrotic gene expression program**, or whether they exhibit **distinct molecular signatures**, is critical for:
- Improving disease classification
- Identifying potential biomarkers
- Informing future therapeutic strategies

---

## Study Design

### Comparisons Performed
Three biologically meaningful contrasts were evaluated:
1. **CHP vs Control** — identifies genes associated with CHP-related fibrosis
2. **IPF vs Control** — identifies genes associated with IPF-related fibrosis
3. **CHP vs IPF** — identifies genes distinguishing the two diseases directly

---

## Methodology (Summary)

### Data Source
- **GEO accession:** GSE150910  
- **Data type:** RNA-seq gene-level raw count matrix  
- **Samples:** Lung tissue from Control, CHP, and IPF patients  

### Data Processing
- Raw count data were used (no pre-normalization)
- Genes with zero counts across samples were removed
- Sample metadata were aligned and validated

### Exploratory Quality Control
- **Variance-stabilizing transformation (VST)** applied
- **PCA** used to assess dominant sources of variation and sample clustering
- **Hierarchical clustering and sample-to-sample distance heatmaps** used to detect outliers and assess biological grouping

### Differential Expression Analysis
- Performed using **DESeq2**
- Negative binomial generalized linear models
- Significance thresholds:
  - False Discovery Rate (FDR) < 0.05
  - |log2 fold change| > 1

### Signature Definition
- **CHP signature:** significant genes in CHP vs Control
- **IPF signature:** significant genes in IPF vs Control
- **Shared fibrotic signature:** genes significant in both contrasts with the same direction of effect
- **Disease-specific genes:** unique to CHP or IPF

---

## Key Findings

### 1. Global Expression Patterns
- PCA and clustering analyses show that **gene expression profiles cluster primarily by disease status**
- Control samples are distinct from fibrotic disease samples
- CHP and IPF show partial overlap, reflecting shared fibrosis biology, but also demonstrate separable molecular features

### 2. Shared Fibrotic Gene Expression Program
- A substantial set of genes is **consistently upregulated in both CHP and IPF**
- These genes include known fibrosis-associated markers involved in:
  - Extracellular matrix remodeling
  - Cell adhesion
  - Tissue scarring
- This supports the existence of a **common fibrotic transcriptional program**

### 3. Disease-Specific Differences
- Direct comparison of CHP vs IPF reveals genes that differ in expression between the two diseases
- These differences suggest **distinct molecular pathways** underlying disease progression despite similar fibrotic outcomes

### 4. Marker Gene Validation
- Boxplots of known fibrosis-related genes (e.g., BGN, COL17A1, CTHRC1) confirm:
  - Low expression in controls
  - Elevated and variable expression in CHP and IPF
- These results are consistent with prior biological knowledge and validate the DE analysis

---

## Outputs Generated

- **PCA plots (VST-transformed counts)**
- **Volcano plots** for all three contrasts
- **Heatmaps** of top differentially expressed genes
- **Sample-to-sample distance heatmaps**
- **Marker gene expression boxplots**
- **Table of shared differentially expressed genes** (analogous to Table 5 in the original paper)

All figures are saved in the `figures/` directory and result tables in `results/`.

---

## Reproducibility
All analysis steps are fully scripted in R. Intermediate objects and final results are saved to ensure that the analysis can be rerun and audited at any stage.

---

## Relationship to the Original Study
This work is a **reproduction and validation study** of the findings reported in:

> *Chronic Hypersensitivity Pneumonitis, an Interstitial Lung Disease with Distinct Molecular Signatures*  
> (American Journal of Respiratory and Critical Care Medicine)

No novel claims are made beyond what can be supported by the reproduced analyses. The intent is to:
- Verify the robustness of the reported molecular signatures
- Demonstrate transparent computational reproducibility
- Provide an educational and methodological reference implementation

All credit for the original biological hypotheses, experimental design, and primary conclusions belongs to the original authors.

---

## Author
**Enock Kumi Ackaah**  
Graduate Researcher — Statistics / Bioinformatics  

---

## License
This repository is provided for **academic and educational purposes only**.

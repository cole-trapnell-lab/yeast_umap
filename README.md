# Dimensionality reduction by UMAP to visualize physical and genetic interactions

## Introduction

Dimensionality reduction is often used to visualize complex expression profiling data. Here, we use the Uniform Manifold Approximation and Projection (UMAP) method on published transcript profiles of 1484 single gene deletions of Saccharomyces cerevisiae. Proximity in low-dimensional UMAP space identifies groups of genes that correspond to protein complexes and pathways, and finds novel protein interactions, even within well-characterized complexes. This approach is more sensitive than previous methods and should be broadly useful as additional transcriptome datasets become available for other organisms.

This is a tutorial to get started working with the Kemmeren et al., yeast deletion strain data. It includes:

1. Dimensionality reduction with UMAP and clustering approaches
2. Differential expression testing 
3. Measuring distance between deletion strains in UMAP space 

### Prerequisites

To get started, you will need to install [Monocle 3] (https://cole-trapnell-lab.github.io/monocle3/)

Other helpful R packages are listed at the top of the notebooks. 


## Authors and Citations

### This study

Michael W. Dorrity, Lauren M. Saunders, Christine Queitsch, Stanley Fields & Cole Trapnell. Dimensionality reduction by UMAP to visualize physical and genetic interactions. Nat Commun 11, 1537 (2020). https://doi.org/10.1038/s41467-020-15351-4

### Source Data

Kemmeren, Patrick, Katrin Sameith, Loes A. L. van de Pasch, Joris J. Benschop, Tineke L. Lenstra, Thanasis Margaritis, Eoghan O’Duibhir, et al. 2014. Large-Scale Genetic Perturbations Reveal Regulatory Networks and an Abundance of Gene-Specific Repressors. Cell 157 (3): 740–52. https://doi.org/10.1016/j.cell.2014.02.054

Welcome to gsMap's documentation!
===================================


Introduction
------------

``gsMap`` (genetically informed spatial mapping of cells for complex traits) integrates spatial transcriptomics (ST) data with genome-wide association study (GWAS) summary statistics to map cells to human complex traits, including diseases, in a spatially resolved manner.


How to Cite
------------
If you use ``gsMap`` in your studies, please cite:

   Liyang Song, Wenhao Chen, Junren Hou, Minmin Guo, Jian Yang (2024) Spatially resolved mapping of cells associated with human complex traits. (Under review)

Key Features
------------

- **Spatially-aware High-Resolution Trait Mapping**: Maps trait-associated cells at single-cell resolution, offering insights into their spatial distributions.
- **Spatial Region Identification**: Aggregates trait-cell association p-values into trait-tissue region association p-values, prioritizing tissue regions relevant to traits of interest.
- **Putative Causal Genes Identification**: Prioritizes putative causal genes by associating gene expression levels with cell-trait relevance.


Overview of ``gsMap`` Method
-----------------------------

``gsMap`` operates on a four-step process:

1. **Gene Specificity Assessment in Spatial Contexts**: To address technical noise and capture spatial correlations of gene expression profiles in ST data, ``gsMap`` leverages GNNs to identify homogeneous spots for each spot and estimates gene specificity scores by aggregating information from those homogeneous spots.
2. **Linking Gene Specificity to SNPs**: ``gsMap`` assigns gene specificity scores to single nucleotide polymorphisms (SNPs) based on their proximity to gene transcription start sites (TSS) and SNP-to-gene epigenetic linking maps.
3. **Spatial S-LDSC**: To estimate the relevance of spots to traits, ``gsMap`` associates stratified LD scores of individual spots with GWAS summary statistics using the S-LDSC framework.
4. **Spatial Region Identification**: To evaluate the association of a specific spatial region with traits, ``gsMap`` employs the Cauchy combination test to aggregate p-values from individual spots within that spatial region.

.. image:: _static/schematic.svg
   :width: 600
   :alt: Model architecture
Schmatics of ``gsMap`` method. For more details about the ``gsMap``, please check out our `publication <URL>`__.

Installation
------------

``gsMap`` is available on `gsMap GitHub <https://github.com/LeonSong1995/gsMap>`__.


How to install ``gsMap``, check out the `installation guide <install.rst>`__

Tutorials
---------
How to use ``gsMap``, check out the `tutorials <tutorials.rst>`__


Online Analysis Service (coming soon)
---------
Users could upload their own GWAS summary statistics data to perform the analysis.


.. toctree::
    :maxdepth: 2
    :caption: Contents:

    install
    tutorials
    data
    api
    release


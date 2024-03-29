Welcome to gsMap's documentation!
===============================


Introduction
------------

``gsMap`` (genetically informed spatial mapping of cells for complex traits) is introduced to integrate spatial transcriptomics (ST) data with summary statistics of genome-wide association studies (GWAS), aiming to identify complex trait (disease)-associated cells and map their spatial distributions.

Key Features
------------

- **GWAS and ST Data Integration**: Combines trait--genetic variant association signals and spatial gene expression data for comprehensive analysis.
- **Addressing Technical Noise and Sparsity**: Deals with technical noise and gene expression sparsity in ST data by using GNN.
- **Spatially-aware High-Resolution Trait Mapping**: Maps trait-associated cells at a spatially resolved single-cell resolution, providing insights into their spatial organizations.
- **Spatial Region Identification**: Aggregates trait-cell association statistics into spatial regions, prioritizing spatial-regions related to traits.

Overview of ``gsMap`` Method
--------------------------

``gsMap`` operates on a four-step process:

1. **Gene Specificity Assessment in Spatial Contexts**: To address technical noise and capture spatial correlations of gene expression profiles in ST data, ``gsMap`` leverages GNN to identify neighboring spots for each spot and estimates the gene specificity score by aggregating information from neighboring spots.

2. **Linking Gene Specificity to SNPs**: ``gsMap`` assigns gene specificity scores to single nucleotide polymorphisms (SNPs), based on their proximity to gene transcription start sites (TSS) and SNP-to-gene epigenetic linking maps. Using SNP annotations from individual spots, ``gsMap`` generates stratified linkage disequilibrium (LD) scores for each spot.

3. **Spatial S-LDSC**: To estimate the relevance of spots to traits, ``gsMap`` associates stratified LD scores of spots with GWAS summary statistics of traits using the framework of S-LDSC.

4. **Spatial Region Identification**: To evaluate the association of a specific spatial region with traits, ``gsMap`` employs the Cauchy combination test to aggregate P values of individual spots within that spatial region.

.. image:: _static/architecture.svg
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

Web Application
---------------

You could visit our `gsMap website <https://gps.yanglab.westlake.edu.cn/>`__ to see the results of our analysis.

Online Analysis Service (not avaliable yet)
+++++++++++++++++++++++

Please check out our `gsMap online application <not available yet>`__ for a user-friendly online analysis. Users could upload their own GWAS summary statistics data to perform the analysis.

How to Cite
-----------

If you use ``gsMap`` in your studies, please cite `[URL] <URL>`__.

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    install
    tutorials
    data
    api
    release


GPS Tutorials
=============

Welcome to the GPS Tutorials. In this section, we present detailed examples and guides to help you understand and utilize GPS effectively. Our tutorials cover two primary datasets: Human Cortex ST data and Mouse Embryo ST data, alongside three cleaned GWAS datasets to demonstrate the application of GPS in analyzing the spatial heterogeneity of complex traits.

Datasets Used in the Tutorials
-------------------------------

1. **Human Cortex ST Data**

   - Num of spots: 4k
   - Technology: 10x Genomics Visium.

2. **Mouse Embryo ST Data**

   - Num of spots: 121k
   - Technology: Stereo-Seq.

3. **GWAS Data**

   - The following Genome-Wide Association Study (GWAS) datasets are provided to illustrate the use of GPS in understanding spatial heterogeneity in complex traits:

     1. Height - (https://www.nature.com/articles/s41586-022-05275-y)
     2. IQ - (https://www.nature.com/articles/s41588-018-0152-6)
     3. Schizophrenia - (https://www.nature.com/articles/s41586-022-04434-5)

For the full dataset used in the paper, refer to the following documentation:

- :ref:`data-availability`



Tutorials
---------------

The GPS workflow is generally consistent across both datasets. However, for the Mouse Embryo ST data, homologs conversion is need in `latent_to_gene` step.

The tutorials for each dataset are provided in separate documents for detailed, step-by-step guidance:

.. toctree::

    cortex.md
    mouse.md

Tutorials
=========

This is the tutorial of GPS. We provide two example ST data: the human cortex and mouse embryo ST data.

Human cortex ST data:
    4k spots, 10x Genomics Visium

Mouse embryo ST data:
    121k spots, Stereo-Seq


We also provide tree cleaned GWAS data illustrate how to use ``GPS`` to show the spatial heterogeneity of the complex traits.

GWAS data:

1. Height https://www.nature.com/articles/s41586-022-05275-y
2. IQ https://www.nature.com/articles/s41588-018-0152-6
3. Schizophrenia https://www.nature.com/articles/s41586-022-04434-5


The workflow of GPS is basically the same for these two datasets, except for the mouse data in the `latent_to_gene` step, we need to use the homologous genes between human and mouse to transfer the

See the :ref:`download-example-data`.

.. toctree::

    cortex.md
    mouse.md


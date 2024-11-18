.. _api-documentation:

gsMap Command Line Alphabets
==============================


.. program:: gsmap

Synopsis
--------

.. code-block:: shell

    usage: gsmap [-h] [--version] {run_find_latent_representations,run_latent_to_gene,run_generate_ldscore,run_spatial_ldsc,run_cauchy_combination,run_report,format_sumstats,quick_mode} ...

Description
-----------

gsmap: genetically informed spatial mapping of cells for complex traits.

Options
-------

.. option:: -h, --help

    Show this help message and exit.

.. option:: --version, -v

    Show program's version number and exit.

Subcommands
-----------

.. option:: format_sumstats

    Convert GWAS summary statistics into the format that gsMap can recognize.

.. option:: run_find_latent_representations

    Find the latent representations of each spot by running GNN-VAE.

.. option:: run_latent_to_gene

    Generate gene specificity scores (GSS) for each spot.

.. option:: run_generate_ldscore

    Generate LD scores for each spot.

.. option:: run_spatial_ldsc

    Perform LDSC for each spot.

.. option:: run_cauchy_combination

    Perform cauchy combination test for each annotation.

.. option:: run_report

    Generate gsMap report.

.. option:: quick_mode

    Run entire gsMap Pipeline. 

.. option:: create_slice_mean

    Create Slice Mean using multiple slices.

-----

.. toctree::
    :maxdepth: 1
    :caption: Subcommands Documentation

    api/format_sumstats
    api/quick_mode
    api/find_latent_representations
    api/latent_to_gene
    api/generate_ldscore
    api/spatial_ldsc
    api/cauchy_combination
    api/report
    api/create_slice_mean


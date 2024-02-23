from pathlib import Path

import numpy as np

workdir: '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/macaque/processed'
# workdir: '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/macaque'
sample_name = "Cortex_151507"
chrom = "all"
QOS = "huge"
# chrom = range(1,23)
trait_names = [
    'PGC3_SCZ_wave3_public_INFO80'
]
root = "/storage/yangjianLab/songliyang/SpatialData/Data/Brain/macaque/Cell/processed/h5ad"
# sample_names = [file.strip().split('.')[0]
#                 for file in open(f'{root}/representative_slices2').readlines()]
#
# sample_names = '''
# T33_macaque1  T44_macaque1  T82_macaque1  T97_macaque1  T125_macaque1  T127_macaque1  T129_macaque1  T131_macaque1  T135_macaque1  T137_macaque1  T139_macaque1
# '''.strip().split()
sample_names=[]
for file in Path(root).glob('*.h5ad'):
    sample_names.append(file.stem)
sample_names.remove('T825_macaque3') # due to 25% of spot don't have spatial coordinates

annotation = "SubClass"
data_type = "SCT"
# sample_names = ['T584_macaque2']

rule all:
    input:
        expand('{sample_name}/spatial_ldsc/{sample_name}.spatial_ldsc.done',trait_name=trait_names,sample_name=sample_names)


# expand('{sample_name}/cauchy_combination/{sample_name}_{trait_name}.Cauchy.csv.gz',trait_name=trait_names,sample_name=sample_names)
    # expand('{sample_name}/cauchy_combination/{sample_name}_{trait_name}.Cauchy.csv.gz',trait_name=trait_names,sample_name=sample_names)

rule test_run:
    input:
        [f'{sample_name}/generate_ldscore/{sample_name}_generate_ldscore_chr{chrom}.done' for sample_name in
         sample_names]

# localrules: find_latent_representations,latent_to_gene
def get_annotation(wildcards):
    if wildcards.sample_name.endswith('3'):
        print(wildcards.sample_name,'will use None as annotation')
        return None
    else:
        print(wildcards.sample_name,'will use SubClass as annotation')
        return 'SubClass'


rule find_latent_representations:
    input:
        hdf5_path=f'{root}/{{sample_name}}.h5ad'
    output:
        hdf5_output='{sample_name}/find_latent_representations/{sample_name}_add_latent.h5ad'
    params:
        annotation= get_annotation,
        type=data_type,
        epochs=300,
        feat_hidden1=256,
        feat_hidden2=128,
        feat_cell=3000,
        gcn_hidden1=64,
        gcn_hidden2=30,
        p_drop=0.1,
        gcn_lr=0.001,
        gcn_decay=0.01,
        n_neighbors=11,
        label_w=1,
        rec_w=1,
        n_comps=300,
        weighted_adj=False,
        nheads=3,
        var=False,
        convergence_threshold=1e-4,
        hierarchically=False
    threads:
        3
    benchmark: '{sample_name}/find_latent_representations/{sample_name}_add_latent.h5ad.benchmark'
    resources:
        mem_mb_per_cpu=lambda wildcards, threads, attempt: 20_000 * np.log2(attempt + 1),
        qos=QOS
    run:
        command = f"""
GPS run_find_latent_representations \
    --input_hdf5_path {input.hdf5_path} \
    --sample_name {wildcards.sample_name} \
    --output_hdf5_path {output.hdf5_output} \
    { '--annotation ' + params.annotation if params.annotation is not None else ''} \
    --type {params.type} \
    --epochs {params.epochs} \
    --feat_hidden1 {params.feat_hidden1} \
    --feat_hidden2 {params.feat_hidden2} \
    --feat_cell {params.feat_cell} \
    --gcn_hidden1 {params.gcn_hidden1} \
    --gcn_hidden2 {params.gcn_hidden2} \
    --p_drop {params.p_drop} \
    --gcn_lr {params.gcn_lr} \
    --gcn_decay {params.gcn_decay} \
    --n_neighbors {params.n_neighbors} \
    --label_w {params.label_w} \
    --rec_w {params.rec_w} \
    --n_comps {params.n_comps} \
    {'--weighted_adj' if params.weighted_adj else ''} \
    --nheads {params.nheads} \
    {'--var' if params.var else ''} \
    --convergence_threshold {params.convergence_threshold} \
    {'--hierarchically' if params.hierarchically else ''}
        """
        shell(
            f'{command}'
        )


rule latent_to_gene:
    input:
        hdf5_with_latent_path=rules.find_latent_representations.output.hdf5_output
    output:
        feather_path='{sample_name}/latent_to_gene/{sample_name}_gene_marker_score.feather'
    params:
        latent_representation="latent_GVAE",
        num_neighbour=51,
        num_neighbour_spatial=201,
        species='MACAQUE_GENE_SYM',
        gs_species='/storage/yangjianLab/songliyang/SpatialData/homologs/macaque_human_homologs.txt',
        gM_slices=None,
        annotation=get_annotation,
        type=data_type
    threads:
        1
    resources:
        mem_mb_per_cpu=lambda wildcards, threads, attempt: 70_000 * np.log2(attempt + 1),
        qos=QOS
    benchmark: '{sample_name}/latent_to_gene/{sample_name}_gene_marker_score.feather.benchmark'
    run:
        command = f"""
GPS run_latent_to_gene \
    --input_hdf5_with_latent_path {input.hdf5_with_latent_path} \
    --sample_name {wildcards.sample_name} \
    --output_feather_path {output.feather_path} \
    { '--annotation ' + params.annotation if params.annotation is not None else ''} \
    --type {params.type} \
    --latent_representation {params.latent_representation} \
    --num_neighbour {params.num_neighbour} \
    --num_neighbour_spatial {params.num_neighbour_spatial} \
     {'--species ' + params.species if params.species is not None else ''} \
     {'--gs_species ' + params.gs_species if params.gs_species is not None else ''} \
     {'--gM_slices ' + params.gM_slices if params.gM_slices is not None else ''}
"""
        shell(
            f'{command}'
        )


rule generate_ldscore:
    input:
        mkscore_feather_file=rules.latent_to_gene.output.feather_path
    output:
        done='{sample_name}/generate_ldscore/{sample_name}_generate_ldscore_chr{chrom}.done'
    params:
        ld_score_save_dir='{sample_name}/generate_ldscore',
        gtf_annotation_file="/storage/yangjianLab/songliyang/ReferenceGenome/GRCh37/gencode.v39lift37.annotation.gtf",
        bfile_root="/storage/yangjianLab/sharedata/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC",
        keep_snp_root="/storage/yangjianLab/sharedata/LDSC_resource/hapmap3_snps/hm",
        gene_window_size=50000,
        enhancer_annotation_file=None,
        snp_multiple_enhancer_strategy='max_mkscore',
        gene_window_enhancer_priority=None,
        spots_per_chunk=1000,
        ld_wind=1,
        ld_unit="CM",
        additional_baseline_annotation_dir_path=None
        # additional_baseline_annotation_dir_path='/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/resource/ldsc/baseline_v1.2/remove_base'
    benchmark: '{sample_name}/generate_ldscore/{sample_name}_generate_ldscore_chr{chrom}.done.benchmark'
    threads:
        3
    resources:
        mem_mb_per_cpu=lambda wildcards, threads, attempt: 50_000 / threads * np.log2(attempt + 1),
        qos=QOS
    run:
        command = f"""
        GPS run_generate_ldscore \
            --sample_name {wildcards.sample_name} \
            --chrom {wildcards.chrom} \
            --ldscore_save_dir {params.ld_score_save_dir} \
            --mkscore_feather_file {input.mkscore_feather_file} \
            --bfile_root {params.bfile_root} \
            --keep_snp_root {params.keep_snp_root} \
            --gtf_annotation_file {params.gtf_annotation_file} \
            --gene_window_size {params.gene_window_size} \
            {'--enhancer_annotation_file ' + params.enhancer_annotation_file if params.enhancer_annotation_file is not None else ''} \
            --snp_multiple_enhancer_strategy {params.snp_multiple_enhancer_strategy} \
            {'--gene_window_enhancer_priority ' + params.gene_window_enhancer_priority if params.gene_window_enhancer_priority is not None else ''} \
            --spots_per_chunk {params.spots_per_chunk} \
            --ld_wind {params.ld_wind} \
            --ld_unit {params.ld_unit} \
            { '--additional_baseline_annotation_dir_path ' +  params.additional_baseline_annotation_dir_path if params.additional_baseline_annotation_dir_path is not None else '' }
        """
        shell(command)
        shell('touch {output.done}')


def get_h2_file(wildcards):
    gwas_root = "/storage/yangjianLab/songliyang/GWAS_trait/LDSC"
    return f"{gwas_root}/{wildcards.trait_name}.sumstats.gz",


def get_ldscore(wildcards):
    if chrom == "all":
        return f"{wildcards.sample_name}/generate_ldscore/{wildcards.sample_name}_generate_ldscore_chr{chrom}.done"
    else:
        assert tuple(chrom) == tuple(range(1,23)), "chrom must be all or range(1,23)"
        return [f"{wildcards.sample_name}/generate_ldscore/{wildcards.sample_name}_generate_ldscore_chr{c}.done" for
                c in chrom]


rule spatial_ldsc:
    input:
        # h2_file=get_h2_file,
        generate_ldscore_done=get_ldscore
    output:
        done='{sample_name}/spatial_ldsc/{sample_name}.spatial_ldsc.done'
    params:
        ldscore_input_dir=rules.generate_ldscore.params.ld_score_save_dir,
        ldsc_save_dir='{sample_name}/spatial_ldsc',
        w_file="/storage/yangjianLab/sharedata/LDSC_resource/LDSC_SEG_ldscores/weights_hm3_no_hla/weights.",
        sumstats_config_file='/storage/yangjianLab/chenwenhao/projects/202312_GPS/src/GPS/example/sumstats_config_sub.yaml',
        all_chunk = None
    threads:
        10
    benchmark:
        '{sample_name}/spatial_ldsc/{sample_name}.spatial_ldsc.done.benchmark'
    resources:
        mem_mb_per_cpu=lambda wildcards, threads, attempt: 40_000 / threads * np.log2(attempt + 1),
        qos=QOS,
        partition='intel-sc3,amd-ep2'
    run:
       command = f"""
        GPS run_spatial_ldsc --w_file {params.w_file} --sample_name {wildcards.sample_name} --num_processes {threads} --ldscore_input_dir {params.ldscore_input_dir} --ldsc_save_dir {params.ldsc_save_dir} --sumstats_config_file {params.sumstats_config_file} {f'--all_chunk {params.all_chunk}' if params.all_chunk else ''}
        """
       shell(
           f'{command}'
           'touch {output.done}'
       )


rule cauchy_combination:
    output:
        done='{sample_name}/cauchy_combination/{sample_name}_{trait_name}.Cauchy.csv.gz'
    input:
        hdf5_path=rules.find_latent_representations.output.hdf5_output,
        ldsc_done=rules.spatial_ldsc.output.done
    params:
        cauchy_save_dir='{sample_name}/cauchy_combination',
        annotation=annotation,
        ldsc_dir=rules.spatial_ldsc.params.ldsc_save_dir
    benchmark:
        '{sample_name}/cauchy_combination/{sample_name}_{trait_name}.Cauchy.csv.gz.benchmark'
    threads:
        2
    resources:
        mem_mb_per_cpu=25_000
    shell:
        """
        GPS run_cauchy_combination --input_hdf5_path {input.hdf5_path} --input_ldsc_dir {params.ldsc_dir} --sample_name {wildcards.sample_name} --output_cauchy_dir {params.cauchy_save_dir} --trait_name {wildcards.trait_name} --annotation {params.annotation}
        """

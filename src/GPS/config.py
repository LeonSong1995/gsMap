from typing import Optional, Literal

from omegaconf import OmegaConf, II
from dataclasses import dataclass, field
from omegaconf import MISSING

@dataclass
class ST_SAMPLE_INFO:
    sample_hdf5:str = MISSING
    sample_name:str = MISSING
    annotation_layer_name: Optional[str] = None
    is_count: bool = True


@dataclass
class FIND_LATENT_REPRESENTATIONS_Config:
    sample_info: ST_SAMPLE_INFO = II("sample_info")
    # epochs: int = 300
    # feat_hidden1: int = 256
    # feat_hidden2: int = 128
    # feat_cell: int = 3000
    # n_comps: int = 300
    # n_neighbors: int = 11
    # p_drop: float = 0.1
    # convergence_threshold: float = 1e-4
    output_dir:str = II('output_dir')


@dataclass
class LATENT_TO_GENE_Config:
    sample_info: ST_SAMPLE_INFO = II("sample_info")
    method: str = 'rank'
    num_neighbour: int = 21
    num_neighbour_spatial: int = 101
    num_processes: int = 4
    fold: float = 1.0
    pst: float = 0.2
    species: Optional[str] = None
    gs_species: Optional[str] = None
    gM_slices: Optional[str] = None
    output_dir:str = II('output_dir')

@dataclass
class MAKE_ANNOTATION_Conifg:
    sample_info: ST_SAMPLE_INFO = II("sample_info")
    gtf_file: str = MISSING
    bfile_root: str = MISSING
    baseline_annotation: Optional[str] = None
    keep_snp_root: Optional[str] = None
    chr: Optional[int] = None
    window_size: int = 50_000
    cells_per_chunk: int = 100
    ld_wind: int = 1
    ld_wind_unit: Literal['CM', 'BP','SNP'] = 'CM'
    r2_cache_dir: Optional[str] = None
    use_gpu: bool = True
    snps_per_chunk :int = 50_000
    output_dir:str = II('output_dir')

@dataclass
class GPS:
    sample_info: ST_SAMPLE_INFO = field(default_factory=ST_SAMPLE_INFO)
    find_latent_representations: FIND_LATENT_REPRESENTATIONS_Config = field(default_factory=FIND_LATENT_REPRESENTATIONS_Config)
    latent_to_gene: LATENT_TO_GENE_Config = field(default_factory=LATENT_TO_GENE_Config)
    make_annotation: MAKE_ANNOTATION_Conifg = field(default_factory=MAKE_ANNOTATION_Conifg)
    output_dir: str = MISSING

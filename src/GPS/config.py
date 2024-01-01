import argparse

from GPS.find_latent_representation import add_find_latent_representations_args, run_find_latent_representation, \
    FindLatentRepresentationsConfig
from GPS.generate_ldscore import add_generate_ldscore_args, run_generate_ldscore, GenerateLDScoreConfig
from GPS.latent_to_gene import add_latent_to_gene_args, run_latent_to_gene, LatentToGeneConfig

parser = argparse.ArgumentParser(description=" GPS: Genetics-informed pathogenic spatial mapping")
def add_find_latent_representations_args(parser):

    parser.add_argument('--epochs', default=300, type=int, help="Number of training epochs for the GNN-VAE model. Default is 300.")
    parser.add_argument('--feat_hidden1', default=256, type=int, help="Number of neurons in the first hidden layer of the feature extraction network. Default is 256.")
    parser.add_argument('--feat_hidden2', default=128, type=int, help="Number of neurons in the second hidden layer of the feature extraction network. Default is 128.")
    parser.add_argument('--feat_cell', default=3000, type=int, help="Number of top variable genes to select. Default is 3000.")
    parser.add_argument('--gcn_hidden1', default=64, type=int, help="Number of units in the first hidden layer of the GCN. Default is 64.")
    parser.add_argument('--gcn_hidden2', default=30, type=int, help="Number of units in the second hidden layer of the GCN. Default is 30.")
    parser.add_argument('--p_drop', default=0.1, type=float, help="Dropout rate used in the GNN-VAE model. Default is 0.1.")
    parser.add_argument('--gcn_lr', default=0.001, type=float, help="Learning rate for the GCN network. Default is 0.001.")
    parser.add_argument('--gcn_decay', default=0.01, type=float, help="Weight decay (L2 penalty) for the GCN network. Default is 0.01.")
    parser.add_argument('--n_neighbors', default=11, type=int, help="Number of neighbors to consider for graph construction in GCN. Default is 11.")
    parser.add_argument('--label_w', default=1, type=float, help="Weight of the label loss in the loss function. Default is 1.")
    parser.add_argument('--rec_w', default=1, type=float, help="Weight of the reconstruction loss in the loss function. Default is 1.")
    parser.add_argument('--input_pca', default=True, type=bool, help="Whether to perform PCA on input features. Default is True.")
    parser.add_argument('--n_comps', default=300, type=int, help="Number of principal components to keep if PCA is performed. Default is 300.")
    parser.add_argument('--weighted_adj', default=False, type=bool, help="Whether to use a weighted adjacency matrix in GCN. Default is False.")
    parser.add_argument('--nheads', default=3, type=int, help="Number of heads in the attention mechanism of the GNN. Default is 3.")
    parser.add_argument('--var', default=False, type=bool)
    parser.add_argument('--convergence_threshold', default=1e-4, type=float, help="Threshold for convergence during training. Training stops if the loss change is below this threshold. Default is 1e-4.")
    parser.add_argument('--hierarchically', default=False, type=bool, help="Whether to find latent representations hierarchically. Default is False.")

    parser.add_argument('--input_hdf5_path', required=True, type=str, help='Path to the input hdf5 file.')
    parser.add_argument('--output_hdf5_path', required=True, type=str, help='Path to the output hdf5 file.')
    parser.add_argument('--sample_name', required=True, type=str, help='Name of the sample.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")

def add_generate_ldscore_args(parser):
    def chrom_choice(value):
        if value.isdigit():
            ivalue = int(value)
            if 1 <= ivalue <= 22:
                return ivalue
        elif value.lower() == 'all':
            return value
        else:
            raise argparse.ArgumentTypeError(f"'{value}' is an invalid chromosome choice. Choose from 1-22 or 'all'.")

    parser.add_argument('--sample_name', type=str, required=True, help='Sample name')
    parser.add_argument('--chrom', type=chrom_choice, required=True, help='Chromosome number (1-22) or "all"')
    parser.add_argument('--save_dir', type=str, required=True, help='Directory to save the data')
    parser.add_argument('--gtf_file', type=str, required=True, help='GTF file path')
    parser.add_argument('--mkscore_feather_file', type=str, required=True, help='Mkscore feather file path')
    parser.add_argument('--bfile_root', type=str, required=True, help='Bfile root path')
    parser.add_argument('--keep_snp_root', type=str, required=True, help='Keep SNP root path')

    # Arguments with defaults
    parser.add_argument('--window_size', type=int, default=50000, help='Annotation window size for each gene')
    parser.add_argument('--spots_per_chunk', type=int, default=10000, help='Number of spots per chunk')
    parser.add_argument('--ld_wind', type=int, default=1, help='LD window size')
    parser.add_argument('--ld_unit', type=str, default='CM', help='LD window unit (SNP/KB/CM)',
                        choices=['SNP', 'KB', 'CM'])

def add_latent_to_gene_args(parser):
    parser.add_argument('--input_hdf5_path', type=str, required=True, help='Path to the input HDF5 file.')
    parser.add_argument('--sample_name', type=str, required=True, help='Name of the sample.')
    parser.add_argument('--output_feather_path', type=str, required=True, help='Path to save output gene marker score feather file.')
    parser.add_argument('--annotation', default=None, type=str, help='Name of the annotation layer.')
    parser.add_argument('--type', default=None, type=str, help="Type of input data (e.g., 'count', 'counts').")

    parser.add_argument('--method', type=str, default='rank', choices=['rank', 'other_method'], help='Method to be used. Default is "rank".')
    parser.add_argument('--latent_representation', type=str, default='latent_GVAE', choices=['latent_GVAE', 'latent_PCA'], help='Type of latent representation. Default is "latent_GVAE".')
    parser.add_argument('--num_neighbour', type=int, default=21, help='Number of neighbours to consider. Default is 21.')
    parser.add_argument('--num_neighbour_spatial', type=int, default=101, help='Number of spatial neighbours to consider. Default is 101.')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of processes to use. Default is 4.')
    parser.add_argument('--fold', type=float, default=1.0, help='Fold change threshold. Default is 1.0.')
    parser.add_argument('--pst', type=float, default=0.2, help='PST value. Default is 0.2.')
    parser.add_argument('--species', type=str, default=None, help='Species name, if applicable.')
    parser.add_argument('--gs_species', type=str, default=None, help='Gene species file path, if applicable.')
    parser.add_argument('--gM_slices', type=str, default=None, help='Path to gene model slices file, if applicable.')

import argparse
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from GPS.config import VisualizeConfig, add_Visualization_args

# The fun of visualization
def visualize(space_coord,plt_pth,title,height,wdith,dpi,text_size,font_size,point_size,fig_facecolor,fig_style):

    fig, ax = plt.subplots(figsize=(height,wdith),dpi=dpi)
    
    if(fig_style=='dark'):
        facecolor = 'black'
        text_color = 'white'
    else:
        facecolor = 'white'
        text_color = 'black'
    
    # Set the background color to black
    fig.set_facecolor(facecolor)
    ax.set_facecolor(fig_facecolor)
    
    # Set the Pvalue color
    custom_colors = ['#313695','#4575b4','#74add1','#abd9e9','#e0f3f8',
                     '#fee090','#fdae61','#f46d43','#d73027','#a50026']
    cmap_custom = mcolors.ListedColormap(custom_colors)

    # Plot the results
    scatter = plt.scatter(space_coord.sx, space_coord.sy, c = space_coord.logp, 
                          cmap = cmap_custom,marker='.',edgecolor='none',lw=0, s=point_size)
    
    # Adjust the cbar
    cbar = plt.colorbar(shrink=0.60,location = 'left')
    cbar.set_label(r"$-\log_{10}({P-value})$",labelpad=2,size=text_size,color=text_color)
    cbar.ax.tick_params(labelsize=text_size, colors=text_color)
    
    # Add labels and title
    # plt.xlabel('Sx')
    # plt.ylabel('Sy')
    plt.title(title,fontsize = font_size,color=text_color)

    # Remove other elements
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Save the plot
    plt.savefig(plt_pth)
    return 1

def run_Visualize(config:VisualizeConfig):
    # Load the ldsc results
    print(f'------Loading LDSC results of {config.input_ldsc_dir}...')
    ldsc_input_file= Path(config.input_ldsc_dir)/f'{config.sample_name}_{config.trait_name}.csv.gz'
    ldsc = pd.read_csv(ldsc_input_file, compression='gzip')
    ldsc.spot = ldsc.spot.astype(str).replace('\.0', '', regex=True)
    ldsc.index = ldsc.spot
    ldsc['logp'] = -np.log10(ldsc.p)

    # Load the ST data
    print(f'------Loading ST data of {config.sample_name}...')
    adata = sc.read_h5ad(f'{config.input_hdf5_path}')

    # Process the ldsc results
    spot_name = adata.obs_names.to_list()
    space_coord = pd.DataFrame(adata.obsm['spatial'],columns=['sx','sy'],index=spot_name)
    space_coord = space_coord[space_coord.index.isin(ldsc.spot)]
    space_coord = pd.concat([space_coord.loc[ldsc.spot],ldsc.logp],axis=1)

    # Visualization
    output_dir = Path(config.output_figure_dir)
    output_dir.mkdir(parents=True, exist_ok=True, mode=0o755)
    output_file = output_dir / f'{config.sample_name}_{config.trait_name}.png'
    
    if config.fig_title is None:
        fig_title = config.trait_name
    else:
        fig_title = config.fig_title

    visualize(space_coord,
              plt_pth=output_file,
              title=fig_title,
              height=config.fig_height,
              wdith=config.fig_width,
              dpi=config.fig_dpi,
              text_size=config.text_size,
              font_size=config.font_size,
              point_size=config.point_size,
              fig_facecolor=config.fig_facecolor,
              fig_style=config.fig_style)
    

if __name__ == '__main__':
    TEST = True
    if TEST:
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'
        name = 'Cortex_151507'

        config = VisualizeConfig(
            input_hdf5_path= f'{test_dir}/{name}/hdf5/{name}_add_latent.h5ad',
            input_ldsc_dir=
            f'/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/snake_workdir/Cortex_151507/ldsc/',
            output_figure_dir='/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/snake_workdir/Cortex_151507/cauchy/',
            sample_name=name,
            trait_name='adult1_adult2_onset_asthma',
            fig_title='Adult Asthma',
            fig_height=6,fig_wdith=7,fig_dpi=300,
            text_size=10,font_size=12,point_size=1,
            fig_facecolor='black',
            fig_style='dark'
        )
    else:
        parser = argparse.ArgumentParser(description="Visualization the results")
        add_Visualization_args(parser)
        args = parser.parse_args()
        config = VisualizeConfig(**vars(args))
        
    run_Visualize(config)
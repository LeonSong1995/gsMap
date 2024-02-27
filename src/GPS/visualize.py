import argparse
from pathlib import Path
from typing import Literal

import scanpy as sc
import numpy as np
import pandas as pd
import plotly.express as px
from GPS.config import VisualizeConfig, add_Visualization_args


def load_ldsc(ldsc_input_file):
    ldsc = pd.read_csv(ldsc_input_file, compression='gzip')
    ldsc.spot = ldsc.spot.astype(str).replace('\.0', '', regex=True)
    ldsc.index = ldsc.spot
    ldsc['logp'] = -np.log10(ldsc.p)
    return ldsc


# %%
def load_st_coord(adata, ldsc, annotation):
    spot_name = adata.obs_names.to_list()
    space_coord = pd.DataFrame(adata.obsm['spatial'], columns=['sx', 'sy'], index=spot_name)
    space_coord = space_coord[space_coord.index.isin(ldsc.spot)]
    space_coord_concat = pd.concat([space_coord.loc[ldsc.spot], ldsc.logp], axis=1)
    space_coord_concat.head()
    annotation = pd.Series(adata.obs[annotation].values, index=adata.obs_names, name='annotation')
    space_coord_concat = pd.concat([space_coord_concat, annotation], axis=1)
    return space_coord_concat


# %%
def draw_scatter(space_coord_concat, title=None, fig_style: Literal['dark', 'light'] = 'light',
                 point_size: int = None, symbol=None, width=800, height=600):
    # change theme to plotly_white
    if fig_style == 'dark':
        px.defaults.template = "plotly_dark"
    else:
        px.defaults.template = "plotly_white"

    custom_color_scale = [
    (1, '#d73027'),  # Red
    (7 / 8, '#f46d43'),  # Red-Orange
    (6 / 8, '#fdae61'),  # Orange
    (5 / 8, '#fee090'),  # Light Orange
    (4 / 8, '#e0f3f8'),  # Light Blue
    (3 / 8, '#abd9e9'),  # Sky Blue
    (2 / 8, '#74add1'),  # Medium Blue
    (1 / 8, '#4575b4'),  # Dark Blue
    (0, '#313695')  # Deep Blue
]
    # custom_color_scale = px.colors.diverging.balance
    custom_color_scale.reverse()
    fig = px.scatter(
        space_coord_concat,
        x='sx',
        y='sy',
        color='logp',  # Color the points by the 'logp' column
        symbol=symbol,
        title=title,
        color_continuous_scale=custom_color_scale,
        range_color=[0, max(space_coord_concat.logp)],
    )

    if point_size is not None:
        fig.update_traces(marker=dict(size=point_size))

    fig.update_layout(legend_title_text='Annotation')

    fig.update_layout(
        autosize=False,
        width=width,
        height=height,
    )

    # Adjusting the legend
    fig.update_layout(
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=-0.39,
            font=dict(
                size=10,
            )
        )
    )
    # change color bar title
    fig.update_layout(coloraxis_colorbar_title='-log10(p)')

    return fig

def run_Visualize(config: VisualizeConfig):
    print(f'------Loading LDSC results of {config.input_ldsc_dir}...')
    ldsc = load_ldsc(ldsc_input_file=Path(config.input_ldsc_dir) / f'{config.sample_name}_{config.trait_name}.csv.gz')

    print(f'------Loading ST data of {config.sample_name}...')
    adata = sc.read_h5ad(f'{config.input_hdf5_path}')

    space_coord_concat = load_st_coord(adata, ldsc, annotation=config.annotation)
    fig = draw_scatter(space_coord_concat,
                       title=config.fig_title,
                       fig_style=config.fig_style,
                       point_size=config.point_size,
                       symbol='annotation',
                       width=config.fig_width,
                       height=config.fig_height,
    )

    # save the figure to html

    # Visualization
    output_dir = Path(config.output_figure_dir)
    output_dir.mkdir(parents=True, exist_ok=True, mode=0o755)
    output_file_html = output_dir / f'{config.sample_name}_{config.trait_name}.html'
    output_file_pdf = output_dir / f'{config.sample_name}_{config.trait_name}.pdf'

    fig.write_html(str(output_file_html))
    fig.write_image(str(output_file_pdf))

    print(f'------The visualization result is saved in a html file: {output_file_html} which can interactively viewed in a web browser and a pdf file: {output_file_pdf}.')

if __name__ == '__main__':
    TEST = True
    if TEST:
        test_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021'
        name = 'E16.5_E1S1'

        config = VisualizeConfig(
            input_hdf5_path = f'/storage/yangjianLab/songliyang/SpatialData/Data/Embryo/Mice/Cell_MOSTA/h5ad/E16.5_E1S1.MOSTA.h5ad',
            input_ldsc_dir =
        f'/storage/yangjianLab/songliyang/SpatialData/Data/Embryo/Mice/Cell_MOSTA/ldsc_enrichment_frac/E16.5_E1S1/',
        output_figure_dir = '/storage/yangjianLab/chenwenhao/projects/202312_GPS/data/GPS_test/Nature_Neuroscience_2021/snake_workdir/Cortex_151507/figure/',
        sample_name = name,
        trait_name = 'GIANT_EUR_Height_2022_Nature',
        fig_title = 'GIANT_EUR_Height_2022_Nature',
        fig_height = 800,
        fig_width = 800,
        fig_style = 'light',
        point_size = 2,
        annotation = 'annotation',
        )
        run_Visualize(config)
    else:
        parser = argparse.ArgumentParser(description="Visualization the results")
        add_Visualization_args(parser)
        args = parser.parse_args()
        config = VisualizeConfig(**vars(args))

        run_Visualize(config)

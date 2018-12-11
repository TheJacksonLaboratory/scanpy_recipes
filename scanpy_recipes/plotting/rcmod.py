from matplotlib import rcParams
from scanpy.settings import set_figure_params
from cmocean import cm as cmo

def update_figure_params():
    set_figure_params(
        scanpy=True,
        dpi=150,
        dpi_save=450,
        frameon=True,
        vector_friendly=True,
        color_map="Reds",
        format='pdf',
        transparent=False,
        ipython_format='png2x'
    )

    # scanpy sets grids on. woof
    rcParams["axes.grid"] = False

import io
import os
import base64
import pkg_resources
from jinja2 import Environment, PackageLoader, select_autoescape
from subprocess import check_call
import matplotlib
import matplotlib.pyplot as plt

import scanpy.api.logging as logg
#from scanpy.plotting.tools.scatterplots import umap
from scanpy.plotting.anndata import scatter

from .utils import datestamp


def fig_to_bytes(fig):
    bytesio = io.BytesIO()
    fig.savefig(bytesio, format="png")
    bytesio.seek(0)
    encoded = base64.b64encode(bytesio.read()).decode("ascii")
    return encoded


def _get_sampleid(adata):
    sampleids = adata.uns.get("sampleids", None)
    return sampleids if sampleids is not None else [adata.uns["sampleid"]]


def _get_output_file(adata, ext="html"):
    outname = f"{adata.uns['sampleid']}_{datestamp()}_report.{ext}"
    return os.path.join(adata.uns["output_dir"], outname)


def _load_resource(filename):
    full_path = pkg_resources.resource_filename("scanpy_recipes", filename)
    contents = ""
    try:
        with open(full_path, "r") as fin:
            contents = fin.read()
    except Exception as e:
        print(e)
    return contents


class SCBLReport(object):
    MIN_PAGE = 1
    MAX_PAGE = 6

    CSS_FILES = ("static/bootstrap.min.css", )
    JS_FILES = ("static/jquery-3.3.1.min.js",
                "static/bootstrap.min.js", )

    def __init__(self, ):
        self.env = Environment(
            loader=PackageLoader("scanpy_recipes", "templates"),
            autoescape=select_autoescape(['html']),
            trim_blocks=True,
            lstrip_blocks=True
        )

    def add_report_figures(
        self,
        adata,
        violins=None,
        scatters=None,
        ranks=None,
        cluster_key="cluster",
        batch_key="sampleid"
    ):

        img_dict = {"qc": dict()}
        sampleids = _get_sampleid(adata)
        adata.uns["sampleids"] = sampleids

        if not violins:
            raise Exception("Need QC violin plot.")
        else:
            if isinstance(violins, matplotlib.figure.Figure):
                violins = [violins]
            assert len(sampleids) == len(violins)
            img_dict["qc"]["violin"] = dict((p1, fig_to_bytes(p2))
                                            for p1, p2 in zip(sampleids, violins))

        if not scatters:
            raise Exception("Need QC scatter plot.")
        else:
            if isinstance(scatters, matplotlib.figure.Figure):
                scatters = [scatters]
            assert len(sampleids) == len(scatters)
            img_dict["qc"]["scatter"] = dict((p1, fig_to_bytes(p2)) for p1, p2 in zip(sampleids, scatters))

        if not ranks:
            raise Exception("Need QC rank plot.")
        else:
            if isinstance(ranks, matplotlib.figure.Figure):
                violins = [ranks]
            assert len(sampleids) == len(ranks)
            img_dict["qc"]["rank"] = dict((p1, fig_to_bytes(p2))
                                           for p1, p2 in zip(sampleids, ranks))

        fig = scatter(adata, basis="umap", color=cluster_key, legend_loc="on data", show=False).figure
        img_dict["clusters"] = fig_to_bytes(fig)
        plt.close()

        if adata.uns.get("is_aggregation", False):
            fig = scatter(adata, basis="umap", color=batch_key, legend_loc="right margin", show=False).figure
            img_dict["batches"] = fig_to_bytes(fig)
            plt.close()

        adata.uns["report_images"] = img_dict


    def _render_page(self, adata, n):
        template = self.env.get_template(f"page{n}.html")
        return template.render(adata=adata)


    def generate_report(self, adata):
        # hack
        del adata.uns["10x_metrics"]["sample"]["Sequencing Saturation"]

        report_template = self.env.get_template("report.html")
        pages = [self._render_page(adata, n)
                 for n in range(self.MIN_PAGE, self.MAX_PAGE + 1)]

        css = "\n".join(_load_resource(css_file) for css_file in self.CSS_FILES)
        js = "\n".join(_load_resource(js_file) for js_file in self.JS_FILES)

        html_report = report_template.render(
            adata=adata,
            css=css,
            js=js,
            html="\n".join(pages)
        )

        report_file = _get_output_file(adata, ext="html")
        with open(report_file, "w") as htmlout:
            htmlout.write(html_report)
        logg.info(f"HTML report saved to [{report_file}].")

        self.html_report = html_report
        self.html_file = report_file


    def generate_pdf(self):
        pdf_file = f"{os.path.splitext(self.html_file)[0]}.pdf"
        cmd = ["pandoc", "-f", "html", "-t", "latex", "+RTS", "-K512m", "-RTS", "-o",
                pdf_file, self.html_file]
        try:
            check_call(cmd)
        except Exception as e:
            logg.error(e)
        else:
            logg.info(f"PDF report saved to [{pdf_file}].")


__api_objects__ = {
    "SCBLReport": SCBLReport,
    "capture_figure": fig_to_bytes,
}

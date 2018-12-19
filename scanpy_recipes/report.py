import io
import os
import base64
from jinja2 import Environment, PackageLoader, select_autoescape
from subprocess import check_call
import matplotlib.pyplot as plt

import scanpy.api.logging as logg

from .utils import datestamp


def fig_to_bytes(fig):
    bytesio = io.BytesIO()
    fig.savefig(bytesio, format="png")
    bytesio.seek(0)
    encoded = base64.b64encode(bytesio.read()).decode("ascii")
    return encoded


def _get_output_file(adata, ext="html"):
    outname = f"{adata.uns['sampleid']}_{datestamp()}_report.{ext}"
    return os.path.join(adata.uns["output_dir"], outname)


class SCBLReport(object):
    MIN_PAGE = 1
    MAX_PAGE = 5
    def __init__(self, ):
        self.env = Environment(
            loader=PackageLoader("scanpy_recipes", "templates"),
            autoescape=select_autoescape(['html']),
            trim_blocks=True,
            lstrip_blocks=True
        )

    def _render_page(self, adata, n):
        template = self.env.get_template(f"page{n}.html")
        return template.render(adata=adata)


    def generate_report(self, adata):
        report_template = self.env.get_template("report.html")
        pages = [self._render_page(adata, n)
                 for n in range(self.MIN_PAGE, self.MAX_PAGE + 1)]

        html_report = report_template.render(
            adata=adata,
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

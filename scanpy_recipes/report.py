import io
import base64
from jinja2 import Environment, PackageLoader, select_autoescape
import matplotlib.pyplot as plt


def fig_to_bytes(fig):
    bytesio = io.BytesIO()
    plt.savefig(bytesio, format="png")
    bytesio.seek(0)
    encoded = base64.b64encode(bytesio.read())
    return encoded


class SCBLReport(object):
    MIN_PAGE = 1
    MAX_PAGE = 5
    def __init__(self, ):
        self.env = Environment(
            PackageLoader("scanpy_recipes", "templates"),
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

        report = report_template.render(html="\n".join(pages))
        return report

__api_objects__ = {
    "SCBLReport": SCBLReport,
    "capture_figure": fig_to_bytes,
}

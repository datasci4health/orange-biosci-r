import numpy as np
import pandas as pd
from importlib.resources import files

from AnyQt.QtWidgets import QComboBox

from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets.settings import Setting
from Orange.widgets import gui
from Orange.data import Table, Domain, ContinuousVariable, StringVariable

import rpy2.robjects as ro
from rpy2.robjects.packages import importr

from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import default_converter




class OWLimmaDifferentialExpression(OWWidget):
    name = "Differential Expression (limma)"
    description = "Compute differential expression using limma via rpy2"
    icon = str(files("orange3bioscir") / "icons/DifferentialExpression.svg")

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        results = Output("Results", Table)

    annotation_type = Setting("")
    group1 = Setting("")
    group2 = Setting("")

    want_main_area = False

    def __init__(self):
        super().__init__()

        self.data = None

        self.annotation_box = gui.widgetBox(self.controlArea, "Grouping")
        self.annotation_combo = QComboBox()
        self.annotation_combo.currentTextChanged.connect(self.annotation_changed)
        self.annotation_box.layout().addWidget(self.annotation_combo)

        self.control_box = gui.widgetBox(self.controlArea, "Control")
        self.group1_combo = QComboBox()
        self.control_box.layout().addWidget(self.group1_combo)

        self.case_box = gui.widgetBox(self.controlArea, "Case")
        self.group2_combo = QComboBox()
        self.case_box.layout().addWidget(self.group2_combo)

        gui.button(self.controlArea, self, "Run", callback=self.run_analysis)

        # Set fixed size for the setup window
        self.setFixedSize(300, 350)

    @Inputs.data
    def set_data(self, data):
        self.data = data
        self.populate_annotations()

    def populate_annotations(self):
        self.annotation_combo.clear()

        if self.data is None:
            return

        annotations = set()

        for var in self.data.domain.attributes:
            for k in var.attributes:
                annotations.add(k)

        self.annotation_combo.addItems(sorted(annotations))

    def annotation_changed(self):

        annotation = self.annotation_combo.currentText()

        values = set()

        for var in self.data.domain.attributes:
            if annotation in var.attributes:
                values.add(var.attributes[annotation])

        self.group1_combo.clear()
        self.group2_combo.clear()

        values = sorted(values)

        self.group1_combo.addItems(values)
        self.group2_combo.addItems(values)

    def run_analysis(self):
        self.error()
        if self.data is None:
            return

        try:
            limma = importr("limma")
            base = importr("base")
            stats = importr("stats")
        except Exception as e:
            self.error(f"Failed to load R packages (limma): {e}")
            return

        annotation = self.annotation_combo.currentText()
        g1 = self.group1_combo.currentText()
        g2 = self.group2_combo.currentText()

        X = self.data.X
        sample_names = [v.name for v in self.data.domain.attributes]

        gene_names = None
        if self.data.domain.metas:
            for meta in self.data.domain.metas:
                if meta.is_string:
                    gene_names = [str(x) for x in self.data.get_column_view(meta)[0]]
                    break
                    
        if gene_names is None:
            gene_names = [f"Gene_{i+1}" for i in range(X.shape[0])]

        samples = []
        groups = []

        for var in self.data.domain.attributes:

            val = var.attributes.get(annotation, None)

            if val == g1 or val == g2:
                samples.append(var.name)
                groups.append(val)

        if len(samples) == 0:
            return

        df = pd.DataFrame(np.array(X), columns=sample_names, index=gene_names)

        df = df[samples]

        with localconverter(default_converter + pandas2ri.converter):
            r_expr = pandas2ri.py2rpy(df)

        ro.globalenv["expr"] = r_expr
        ro.globalenv["groups"] = ro.StrVector(groups)
        ro.globalenv["g1_level"] = g1
        ro.globalenv["g2_level"] = g2

        ro.r("""
        groups_factor <- factor(groups, levels=c(g1_level, g2_level))
        design <- model.matrix(~ groups_factor)
        fit <- lmFit(expr, design)
        fit <- eBayes(fit)
        res <- topTable(fit, coef=2, number=Inf, adjust.method="BH")
        """)

        with localconverter(default_converter + pandas2ri.converter):
            res = ro.r("res")

        res = res.reset_index()

        logfc = res["logFC"].values
        pval = res["P.Value"].values
        adjp = res["adj.P.Val"].values

        domain = Domain([
            ContinuousVariable("logFC"),
            ContinuousVariable("p_value"),
            ContinuousVariable("adj_p_value")
        ], metas=[StringVariable("gene")])

        metas = np.array(res["index"].values, dtype=object).reshape(-1,1)

        X = np.vstack([logfc, pval, adjp]).T

        table = Table(domain, X, metas=metas)

        self.Outputs.results.send(table)

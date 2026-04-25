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

limma = importr("limma")
base = importr("base")
stats = importr("stats")


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

        box = gui.widgetBox(self.controlArea, "Grouping")

        self.annotation_combo = QComboBox()
        self.annotation_combo.currentTextChanged.connect(self.annotation_changed)
        box.layout().addWidget(self.annotation_combo)

        self.group1_combo = QComboBox()
        self.group2_combo = QComboBox()

        box.layout().addWidget(self.group1_combo)
        box.layout().addWidget(self.group2_combo)

        gui.button(self.controlArea, self, "Run", callback=self.run_analysis)

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

        if self.data is None:
            return

        annotation = self.annotation_combo.currentText()
        g1 = self.group1_combo.currentText()
        g2 = self.group2_combo.currentText()

        X = self.data.X
        genes = [v.name for v in self.data.domain.attributes]

        samples = []
        groups = []

        for var in self.data.domain.attributes:

            val = var.attributes.get(annotation, None)

            if val == g1 or val == g2:
                samples.append(var.name)
                groups.append(val)

        if len(samples) == 0:
            return

        df = pd.DataFrame(X.T, columns=genes)

        df = df[samples]

        group_vector = pd.Categorical(groups)

        with localconverter(default_converter + pandas2ri.converter):
            r_expr = pandas2ri.py2rpy(df)
        r_groups = ro.FactorVector(groups)

        ro.globalenv["expr"] = r_expr
        ro.globalenv["groups"] = r_groups

        ro.r("""
        design <- model.matrix(~ groups)
        fit <- lmFit(expr, design)
        fit <- eBayes(fit)
        res <- topTable(fit, coef=2, number=Inf, adjust.method="BH")
        """)

        with localconverter(default_converter + pandas2ri.converter):
            res = pandas2ri.rpy2py(ro.r("res"))

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

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
    auto_commit = Setting(False)

    want_main_area = False

    def __init__(self):
        super().__init__()

        self.data = None

        self.annotation_box = gui.widgetBox(self.controlArea, "Grouping")
        self.annotation_combo = QComboBox()
        self.annotation_combo.currentTextChanged.connect(self.on_annotation_changed)
        self.annotation_box.layout().addWidget(self.annotation_combo)

        self.control_box = gui.widgetBox(self.controlArea, "Control")
        self.group1_combo = QComboBox()
        self.control_box.layout().addWidget(self.group1_combo)

        self.case_box = gui.widgetBox(self.controlArea, "Case")
        self.group2_combo = QComboBox()
        self.group2_combo.currentTextChanged.connect(self.on_group2_changed)
        self.case_box.layout().addWidget(self.group2_combo)
        
        self.group1_combo.currentTextChanged.connect(self.on_group1_changed)

        run_box = gui.hBox(self.controlArea)
        gui.button(run_box, self, "Run", callback=self.commit)
        gui.checkBox(run_box, self, "auto_commit", "Run automatically", callback=self.trigger_commit)

        # Set fixed size for the setup window
        self.setFixedSize(300, 350)

    @Inputs.data
    def set_data(self, data):
        self.data = data
        self.populate_annotations()
        self.trigger_commit()

    def populate_annotations(self):
        self.annotation_combo.blockSignals(True)
        self.annotation_combo.clear()

        if self.data is None:
            self.annotation_combo.blockSignals(False)
            return

        annotations = set()
        for var in self.data.domain.attributes:
            for k in var.attributes:
                annotations.add(k)

        items = sorted(annotations)
        self.annotation_combo.addItems(items)
        if self.annotation_type in items:
            self.annotation_combo.setCurrentText(self.annotation_type)
        else:
            if items:
                self.annotation_type = items[0]
                self.annotation_combo.setCurrentText(self.annotation_type)
            else:
                self.annotation_type = ""
        self.annotation_combo.blockSignals(False)
        self.annotation_changed()

    def on_annotation_changed(self, text):
        if self.annotation_type != text:
            self.annotation_type = text
            self.annotation_changed()
            self.trigger_commit()

    def annotation_changed(self):
        if self.data is None:
            return

        annotation = self.annotation_type

        values = set()
        for var in self.data.domain.attributes:
            if annotation in var.attributes:
                values.add(var.attributes[annotation])

        items = sorted(values)

        self.group1_combo.blockSignals(True)
        self.group1_combo.clear()
        self.group1_combo.addItems(items)
        if self.group1 in items:
            self.group1_combo.setCurrentText(self.group1)
        elif items:
            self.group1 = items[0]
            self.group1_combo.setCurrentText(self.group1)
        self.group1_combo.blockSignals(False)

        self.group2_combo.blockSignals(True)
        self.group2_combo.clear()
        self.group2_combo.addItems(items)
        if self.group2 in items:
            self.group2_combo.setCurrentText(self.group2)
        elif items:
            self.group2 = items[1] if len(items) > 1 else items[0]
            self.group2_combo.setCurrentText(self.group2)
        self.group2_combo.blockSignals(False)

    def on_group1_changed(self, text):
        if self.group1 != text:
            self.group1 = text
            self.trigger_commit()

    def on_group2_changed(self, text):
        if self.group2 != text:
            self.group2 = text
            self.trigger_commit()

    def trigger_commit(self):
        if getattr(self, "auto_commit", False):
            self.commit()

    def commit(self):
        self.error()
        if self.data is None:
            return

        try:
            limma = importr("limma")
            base = importr("base")
            stats = importr("stats")
        except Exception:
            try:
                import os
                os.environ["R_LIBS_USER"] = os.path.expanduser("~/R/library")
                ro.r(f'.libPaths(c("{os.path.expanduser("~/R/library")}", .libPaths()))')
                limma = importr("limma")
                base = importr("base")
                stats = importr("stats")
            except Exception as e:
                self.error(f"Failed to load R packages (limma): {e}. Please ensure limma is installed and try setting the R_LIBS_USER environment variable.")
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
        res <- topTable(fit, coef=2, number=Inf, adjust.method="BH", sort.by="none")
        """)

        with localconverter(default_converter + pandas2ri.converter):
            res = ro.r("res")

        logfc = res["logFC"].values
        pval = res["P.Value"].values
        adjp = res["adj.P.Val"].values

        new_vars = [
            ContinuousVariable("logFC"),
            ContinuousVariable("p_value"),
            ContinuousVariable("adj_p_value")
        ]

        new_domain = Domain(self.data.domain.attributes,
                            self.data.domain.class_vars,
                            self.data.domain.metas + tuple(new_vars))

        new_metas_data = np.vstack([logfc, pval, adjp]).T
        
        if self.data.metas.shape[1] > 0:
            metas_new = np.hstack([self.data.metas, new_metas_data])
        else:
            metas_new = new_metas_data

        table = Table(new_domain, self.data.X, self.data.Y, metas_new)

        self.Outputs.results.send(table)

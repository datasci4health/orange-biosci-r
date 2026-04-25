from importlib.resources import files

NAME = "BioSciR"
DESCRIPTION = "Collection of Bio Science processing widgets interacting with R."
ICON = str(files("orange3bioscir") / "icons" / "BioSciR-category.svg")
BACKGROUND = "#9FFFBD"

PRIORITY = 3

WIDGETS = [
    'OWLimmaDifferentialExpression'
]

# The .py file where each widget is implemented
WIDGET_HELP_PATH = (
    # You can link to documentation here
)

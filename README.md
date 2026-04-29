# orange-biosci-r
Orange widgets for bioinformatics processing using R packages.

# Installation Guide for Differential Expression Widget (limma + R)

This guide explains how to install R and required packages.

## Requirements
- R ≥ 4.2
- R packages: BiocManager, limma
- Python: rpy2

---

## Linux (Ubuntu/Debian)

### Install R
sudo apt update
sudo apt install r-base r-base-dev

### User library
mkdir -p ~/R/library
echo '.libPaths("~/R/library")' >> ~/.Rprofile

### Install packages
R -e ".libPaths('~/R/library'); install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install('limma', ask=FALSE, update=FALSE)"

### Environment variable
echo 'export R_LIBS_USER=~/R/library' >> ~/.bashrc
source ~/.bashrc

### Test
R -e ".libPaths('~/R/library'); library(limma); cat('OK\n')"

---

## Windows

### Install R
https://cran.r-project.org/

### Install packages (R console)
install.packages("BiocManager")
BiocManager::install("limma", ask = FALSE, update = FALSE)

### Test
library(limma)

---

## macOS

### Install R
https://cran.r-project.org/

### User library
mkdir -p ~/R/library
echo '.libPaths("~/R/library")' >> ~/.Rprofile

### Install packages
R -e ".libPaths('~/R/library'); install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install('limma', ask=FALSE, update=FALSE)"

### Env
echo 'export R_LIBS_USER=~/R/library' >> ~/.zshrc
source ~/.zshrc

### Test
R -e ".libPaths('~/R/library'); library(limma); cat('OK\n')"

---

## Python
pip install rpy2

If issues:
pip install --no-binary rpy2 rpy2

---

## Final check
R -e "library(limma)"
python -c "import rpy2.robjects as ro; print('OK')"

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
~~~bash
sudo apt update
sudo apt install r-base r-base-dev
~~~

### User library
~~~bash
mkdir -p ~/R/library
echo '.libPaths("~/R/library")' >> ~/.Rprofile
~~~

### Install packages
~~~bash
R -e ".libPaths('~/R/library'); install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install('limma', ask=FALSE, update=FALSE)"
~~~

### Environment variable
~~~bash
echo 'export R_LIBS_USER=~/R/library' >> ~/.bashrc
source ~/.bashrc
~~~

### Test
~~~bash
R -e ".libPaths('~/R/library'); library(limma); cat('OK\n')"
~~~

---

## Windows

### Install R
https://cran.r-project.org/

### Install packages (R console)
~~~bash
install.packages("BiocManager")
BiocManager::install("limma", ask = FALSE, update = FALSE)
~~~

### Test
~~~bash
library(limma)
~~~

---

## macOS


> This guide covers installing R, BiocManager, and limma on macOS.

---

## 1. Install R



1. Go to [https://cran.r-project.org/](https://cran.r-project.org/)
2. Click **"Download R for macOS"**
3. Download the latest `.pkg` installer (choose the correct version for your chip):
   - **Apple Silicon (M1/M2/M3/M4):** `R-x.x.x-arm64.pkg`
   - **Intel Mac:** `R-x.x.x-x86_64.pkg`
4. Open the `.pkg` file and follow the installation wizard
5. Verify by opening **Terminal** and running:

```bash
R --version
```

---

## 2. Install BiocManager

Open **R (run R.app)** or **RStudio** and run the following in the console:

```r
# Install BiocManager from CRAN
install.packages("BiocManager")

# Verify installation
library(BiocManager)
BiocManager::version()
```

---

## 3. Install limma

```r
# Install limma from Bioconductor
BiocManager::install("limma")

# Verify installation
library(limma)
packageVersion("limma")
```

---


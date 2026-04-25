# Widget Development Environment

## Setting up the Environment

### Creating a Virtual Environment

~~~bash
python3 -m venv .venv
~~~

### Running a Virtual Environment

~~~bash
source .venv/bin/activate
~~~

### Installing Orange for Development

~~~bash
pip install orange3
pip install pyqt5
~~~

#### R environment
~~~bash
sudo apt install r-base

sudo apt update
sudo apt install --no-install-recommends software-properties-common dirmngr
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys '51716619E084DAB9'
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/"

mkdir ~/R
mkdir ~/R/library

R -e ".libPaths('~/R/library'); install.packages('BiocManager', repos='https://cloud.r-project.org/'); BiocManager::install('limma', ask=FALSE, update=FALSE, force=TRUE)"

pip install rpy2

export R_LIBS_USER=~/R/library
~~~

~~~bash
pip install -e .
~~~

~~~bash
orange-canvas
~~~
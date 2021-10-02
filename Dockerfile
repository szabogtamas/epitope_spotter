FROM szabogtamas/jupy_rocker

RUN sudo apt-get update -y && \
    sudo apt-get install -y libxt-dev && \
    sudo apt-get install -y libx11-dev && \
    sudo apt-get install -y libcairo2-dev && \
    sudo apt-get install -y libxml2-dev && \
    sudo apt-get install -y libbz2-dev && \
    sudo apt-get install -y liblzma-dev && \
    sudo apt-get install -y default-jdk && \
    sudo apt-get install -y libglpk-dev

ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c bioconda mhcflurry
RUN conda install -c conda-forge matplotlib numpy pandas seaborn

RUN conda create -n parasail python=3.7 libgcc-ng=9.3
RUN conda install -c bioconda -n parasail parasail-python

RUN conda create -n fred python=2.7
RUN conda install -c bioconda -n fred fred2

ADD ./third_party /usr/local/lib/third_party

RUN mkdir -p /usr/cbs/packages && \
  mkdir -p /usr/cbs/packages/netMHC/4.0 && \
  tar -xvzf /usr/local/lib/third_party/netMHC-4.0a.Linux.tar.gz -C /usr/cbs/packages/netMHC/4.0 && \
  mkdir -p /usr/cbs/bio/src/netMHCII-2.3 && \
  tar -xvzf /usr/local/lib/third_party/netMHCII-2.3.Linux.tar.gz -C /usr/cbs/bio/src/netMHCII-2.3 && \
  mkdir -p /usr/cbs/packages/netMHCcons/1.1/netMHCcons-1.1 && \
  tar -xvzf /usr/local/lib/third_party/netMHCcons-1.1a.Linux.tar.gz -C /usr/cbs/packages/netMHCcons/1.1/netMHCcons-1.1 && \
  mkdir -p /usr/cbs/packages/netchop/3.1 && \
  tar -xvzf /usr/local/lib/third_party/netchop-3.1d.Linux.tar.gz -C /usr/cbs/packages/netchop/3.1

ENV PATH=/usr/local/bin:$PATH

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    optparse \
    docstring \
    plotly \
    heatmaply \
    RColorBrewer \
    ggsci \
    rJava \
    ggridges \
    openxlsx \
    readxl

RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('msa')"

RUN R -e "devtools::install_github('kassambara/ggpubr')"
RUN R -e "devtools::install_github('masato-ogishi/Repitope')"

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json
FROM szabogtamas/jupy_rocker

ARG HOME=""

RUN sudo apt-get update -y && \
    sudo apt-get install -y libxt-dev && \
    sudo apt-get install -y libx11-dev && \
    sudo apt-get install -y libcairo2-dev && \
    sudo apt-get install -y libxml2-dev && \
    sudo apt-get install -y libbz2-dev && \
    sudo apt-get install -y liblzma-dev && \
    sudo apt-get install -y default-jdk && \
    sudo apt-get install -y libglpk-dev && \
    sudo apt-get install -y tcsh

ENV CONDA_DIR $HOME/miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c conda-forge jupyterlab jupytext
RUN conda install -c conda-forge matplotlib numpy pandas seaborn
RUN conda install -c conda-forge dash dash_cytoscape
RUN conda install -c plotly plotly jupyter-dash
RUN conda install -c bioconda mhcflurry
RUN conda install -c bioconda blast

RUN conda create -n parasail python=3.7 libgcc-ng=9.3
RUN conda install -c bioconda -n parasail parasail-python

RUN conda create -n fred python=2.7
RUN conda install -c bioconda -n fred fred2

ADD ./third_party /usr/local/lib/third_party/CBS
RUN mkdir -p /usr/cbs/packages && \
  tar -xvzf /usr/local/lib/third_party/CBS/netMHC-4.0a.Linux.tar.gz -C /usr/cbs/packages && \
  sed -i 's#/usr/cbs/packages/netMHC/4.0/netMHC-4.0#/usr/cbs/packages/netMHC-4.0#g' /usr/cbs/packages/netMHC-4.0/netMHC && \
  sudo chmod -R 777 /usr/cbs/packages/netMHC-4.0/
RUN mkdir -p /usr/cbs/packages && \
  tar -xvzf /usr/local/lib/third_party/CBS/netMHCII-2.3.Linux.tar.gz -C /usr/cbs/packages && \
  sed -i 's#/usr/cbs/bio/src/netMHCII-2.3#/usr/cbs/packages/netMHCII-2.3#g' /usr/cbs/packages/netMHCII-2.3/netMHCII-2.3 && \
  sudo chmod -R 777 /usr/cbs/packages/netMHCII-2.3/
RUN mkdir -p /scratch && \
    chmod 1777 /scratch

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json

# Update default Jupyter launch to use conda base env
RUN echo '#!/bin/bash \
      \n cd /home/rstudio \
      \n source /miniconda3/etc/profile.d/conda.sh \
      \n conda activate base \
      \n jupyter lab --ip=0.0.0.0 --port=8989 --allow-root' \
      > /etc/services.d/jupyter/run
FROM szabogtamas/jupy_rocker

RUN sudo apt-get update -y && \
    sudo apt-get install -y libxt-dev && \
    sudo apt-get install -y libx11-dev && \
    sudo apt-get install -y libcairo2-dev && \
    sudo apt-get install -y libxml2-dev && \
    sudo apt-get install -y libbz2-dev && \
    sudo apt-get install -y liblzma-dev

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    mkdir /root/.conda && \
    bash Miniconda3-latest-Linux-x86_64.sh -b && \
    rm -f Miniconda3-latest-Linux-x86_64.sh

RUN pip3 install numpy && \
    pip3 install pandas && \
    pip3 install matplotlib && \
    pip3 install seaborn

RUN pip3 install mhcflurry

ADD ./third_party /usr/local/lib/third_party

RUN cd /usr/local/lib/third_party && \
  tar -xvzf netchop-3.1d.Linux.tar.gz && \
  /usr/cbs/packages/netchop/3.1/netchop-3.1 && \
  mv netchop-3.1 /usr/cbs/packages/netchop/3.1/netchop-3.1 && \
  cd /home/rstudio/

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
    ggridges \
    pROC \
    openxlsx \
    readxl

RUN R -e "devtools::install_github('kassambara/ggpubr')"

RUN chmod a+rwx -R /home/rstudio

ADD ./configs/rstudio-prefs.json /home/rstudio/.config/rstudio/rstudio-prefs.json
FROM wwliao/netmhc

RUN apt-get update && apt-get install -y tcsh gawk
RUN wget https://downloads.iedb.org/tools/mhci/3.1.1/IEDB_MHC_I-3.1.1.tar.gz
RUN tar -zxvf IEDB_MHC_I-3.1.1.tar.gz

RUN apt-get update && apt-get install -y tcsh gawk
RUN wget https://downloads.iedb.org/tools/mhcii/3.1.2/IEDB_MHC_II-3.1.2.tar.gz
RUN tar -zxvf IEDB_MHC_II-3.1.2.tar.gz

RUN pip install mhcflurry

RUN pip3 install numpy && \
    pip3 install pandas && \
    pip3 install matplotlib && \
    pip3 install seaborn && \
    pip3 install awscli && \
    pip3 install pyBigWig

RUN conda install -c bioconda parasail-python
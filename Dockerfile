FROM docker pull wwliao/netmhc

RUN apt-get update && apt-get install -y tcsh gawk
RUN wget https://downloads.iedb.org/tools/mhci/3.1.1/IEDB_MHC_I-3.1.1.tar.gz
RUN tar -zxvf IEDB_MHC_I-3.1.1.tar.gz

RUN apt-get update && apt-get install -y tcsh gawk
RUN wget https://downloads.iedb.org/tools/mhcii/3.1.2/IEDB_MHC_II-3.1.2.tar.gz
RUN tar -zxvf IEDB_MHC_II-3.1.2.tar.gz

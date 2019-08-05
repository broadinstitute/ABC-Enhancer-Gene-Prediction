FROM ubuntu:18.04

RUN apt-get update
RUN apt-get install -y python python3 virtualenv python-pip python3-pip

WORKDIR /usr/src/app

# Setup the python requirements
COPY requirements2.txt ./
COPY requirements3.txt ./

RUN apt-get install -y zlib1g-dev zlib1g
RUN apt-get install -y  libbz2-dev  
RUN apt-get install -y liblzma-dev

RUN pip2 install --no-cache-dir -r requirements2.txt
RUN pip3 install --no-cache-dir -r requirements3.txt

# Setup samtools
RUN apt-get install -y wget
RUN apt-get install -y libncurses5-dev
RUN wget -O samtools-0.1.19.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download
RUN tar xjf samtools-0.1.19.tar.bz2
WORKDIR /usr/src/app/samtools-0.1.19
RUN make -j 12
ENV PATH=/usr/src/app/samtools-0.1.19/:${PATH}
WORKDIR /usr/src/app/

# Setup bedtools
RUN wget -O bedtools-2.26.0.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
RUN tar xzf bedtools-2.26.0.tar.gz
WORKDIR /usr/src/app/bedtools2/
RUN make -j 12
ENV PATH=/usr/src/app/bedtools2/bin/:${PATH}
WORKDIR /usr/src/app/

# Install MACS2
RUN pip install MACS2 
#RUN git clone https://github.com/taoliu/MACS.git && git checkout macs2python3


# Copy the required scripts
COPY src/ src/

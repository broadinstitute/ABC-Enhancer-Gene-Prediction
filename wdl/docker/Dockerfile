FROM ubuntu:18.04

# Install required libraries
RUN apt-get update && apt-get install -y python python3 virtualenv python-pip python3-pip  zlib1g-dev zlib1g libbz2-dev liblzma-dev wget libncurses5-dev

# Set the working directory
WORKDIR /usr/src/app

# Setup the python requirements
COPY requirements2.txt ./
COPY requirements3.txt ./
RUN pip2 install --no-cache-dir -r requirements2.txt
RUN pip3 install --no-cache-dir -r requirements3.txt

# Setup samtools
RUN wget -O samtools-0.1.19.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download &&  tar xjf samtools-0.1.19.tar.bz2 && cd  /usr/src/app/samtools-0.1.19 &&  make -j 4

# Update the path
ENV PATH=/usr/src/app/samtools-0.1.19/:${PATH}

# Setup bedtools
RUN wget -O bedtools-2.26.0.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && tar xzf bedtools-2.26.0.tar.gz && cd /usr/src/app/bedtools2/ &&  make -j 12

# Update the path
ENV PATH=/usr/src/app/bedtools2/bin/:${PATH}

# Install python packages
RUN pip2 install MACS2 && pip2 install progressbar &&  pip3 install progressbar

# Copy the required scripts
COPY src/ src/

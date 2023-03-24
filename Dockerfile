FROM ubuntu:20.04

RUN apt update

# need these packages to download and build samtools:
# https://github.com/samtools/samtools/blob/1.9/INSTALL
RUN apt install -y wget gcc libz-dev ncurses-dev libbz2-dev liblzma-dev \
    libcurl3-dev libcrypto++-dev make
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar jxf samtools-1.9.tar.bz2 && \
    cd samtools-1.9 && ./configure && make install

CMD ["samtools"]

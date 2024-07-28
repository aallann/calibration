# OS Ubuntu base image
FROM ubuntu:22.04

LABEL maintainer="ostos" email="ucabjo0@ucl.ac.uk"

# suppress interactive CLI prompts
ENV DEBIAN_FRONTEND=noninteractive

# Dependencies 
RUN apt update && apt install -y \
    build-essential \
    git \
    cmake \
    libeigen3-dev \
    libboost-all-dev \
    valgrind \
    liblapack-dev \
    liblapacke-dev \
    libopenblas-dev \
    liblua5.3-dev \
    python3 \
    python3-pip \
    && apt clean \
    && rm -rf /var/lib/apt/lists/*

# dir for bash scripts
RUN mkdir -p /usr/local/bin

# copy scripts
COPY ./scripts/scan.sh /usr/local/bin
RUN chmod +x /usr/local/bin/scan.sh
RUN /usr/local/bin/scan.sh

WORKDIR /usr/calibration
COPY . /usr/calibration/

RUN cmake -S . -B build

RUN cmake --build build -t preprocessing

CMD ["./bin/preprocessing"]
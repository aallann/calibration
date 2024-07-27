FROM ubuntu:22.04


LABEL maintainer = "ostos" email = "ucabjo0@ucl.ac.uk"

# dependencies
RUN apt-get update 

RUN apt-get install -y  build-essential \
    cmake \
    git \
    valgrind \
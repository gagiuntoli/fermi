FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
 git \
 build-essential \
 cmake \
 clang-format

COPY . /fermi

WORKDIR /fermi


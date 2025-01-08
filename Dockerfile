FROM ubuntu:latest

RUN apt-get update && apt-get install -y \
 git \
 build-essential \
 cmake \
 clang-format \
 nodejs \
 npm \
 && rm -rf /var/lib/apt/lists/*

COPY . /fermi

WORKDIR /fermi


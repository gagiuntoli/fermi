FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
 git \
 build-essential \
 cmake \
 clang-format \
 nodejs \
 npm \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /fermi


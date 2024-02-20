#!/bin/bash

clang-format -i $(find src include -name '*.c' -o -name '*.h') 

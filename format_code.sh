#!/bin/bash

clang-format -i $(find src include -name '*.cpp' -o -name '*.hpp') 

#!/usr/bin/env python2

import re
import sys

regex1 = r"(c)(\d+)(-\w+)"
regex2 = r"({)(\w+)"

f = open('channel_name.sh', 'w')
f.write('channel_name=( ')

with open(sys.argv[1], 'r') as test_cases:
    for test in test_cases:
        if re.search(regex1, test):
            match  = re.search(regex1, test)            
            f.write(match.group(0))
            f.write(' ')

f.write(')')

f = open('surface_ID.sh', 'w')
f.write('surface_ID=( ')

with open(sys.argv[1], 'r') as test_cases:
    for test in test_cases:
        if re.search(regex1, test):
            match  = re.search(regex2, test)            
            f.write(match.group(2))
            f.write(' ')

f.write(')')
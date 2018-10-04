#!/usr/bin/python
import re

data = []
i = 0
for p1 in range(10, 110, 10):
    for p2 in range(10, 110, 10):
        data.append([p1, p2 * p1])

j = 0
file = open("result2.txt")
for line in file:
    i += 1
    if i % 11 == 4:
        print line
        mfloat = re.search('.*([0-9]\.[0-9]+)', line).group(1)
        data[j].append(mfloat)
        j += 1

print data

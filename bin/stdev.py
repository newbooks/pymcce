#!/usr/bin/env python

import sys
import numpy as np

if __name__ == "__main__":

    fn_fort38 = sys.argv[1]
    lines = open(fn_fort38).readlines()
    n = len(lines[0].split())
    if n < 1:
        print("fort.38 has an empty line on top. Please remove it.")
        sys.exit()

    confnames = []
    for line in lines[1:]:
        fields = line.split()
        if len(fields) >= 2:
            confnames.append(fields[0])

    fort38 = np.loadtxt(fn_fort38, skiprows=1, usecols=range(1,n)) # for odd reasons, fort.38 has an

    # find standard deviation in each row
    fort38_mean = fort38.mean(axis=1)
    fort38_stdev = fort38.std(axis=1)

    rows, columns = fort38.shape
    print("Conformer       Mean  Stdev")
    for irow in range(rows):
        print("%s %5.3f  %5.3f" %(confnames[irow], fort38_mean[irow], fort38_stdev[irow]))

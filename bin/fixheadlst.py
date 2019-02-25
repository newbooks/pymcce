#!/usr/bin/env python

import sys
import numpy as np
fn_fort38 = "fort.38"

if __name__ == "__main__":

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
    fort38_mean = fort38.mean(axis=1)

    head3lst = open("head3.lst").readlines()
    print(head3lst[0])
    for i in range(len(head3lst)-1):
        line = head3lst[i+1]
        fields = line.split()
        confname = confnames[i]
        if fields[1] == confname:
            pass
        else:
            print("Mismatch %s and %s in line %d of head3lst" %(fields[0], confname, i+1))
            sys.exit()
        if abs(fort38_mean[i]) < 0.0001 or abs(fort38_mean[i]-1.0) < 0.0001:
            FL = "t"
        else:
            FL = fields[2]
        print("%s%s%s" %(line[:21],FL,line[22:].strip("\n")))

#!/usr/bin/env python
"""
This program reads in mc_out and entropy.out to remove entropy from running energy.
It imports entropy.out and fort.38, multiplies them, and writes out the entropy correction at each pH.
This number can be used to correct the "inflated" running and average energy in mc_out.
"""

import sys
import numpy as np

if __name__ == "__main__":
    line = open("fort.38").readline()
    n = len(line.split())
    if n <1:
        print("fort.38 has an empty line on top. Please remove it.")
        sys.exit()

    fort38 = np.loadtxt("fort.38", skiprows=1, usecols=range(1,n)) # for odd reasons, fort.38 has an
    # empty line
    entropy = np.loadtxt("entropy.out", skiprows = 1, usecols=range(1,n))
    result = fort38 * entropy
    entropy_sum = np.sum(result, axis=0)

#    print fort38
    print(line),
    sum_str = " ".join(["%5.2f" % x for x in entropy_sum])
    print("WeightedEntropy%s" % sum_str)
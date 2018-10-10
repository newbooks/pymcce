#!/usr/bin/env python

import sys
import re
from energy import *

if __name__ == "__main__":
    ph = 0.0
    eh = 0.0

    if len(sys.argv) > 1:
        cutoff = float(sys.argv[1])
    else:
        cutoff = PW_PRINT_CUT

    lines = open("states").readlines()
    for line in lines:
        state = [int(x) for x in re.findall(r"[\w']+", line)]
        analyze_state_energy(state, ph=7.0, eh=0.0, cutoff=cutoff)

#    for conf in head3lst:
#        conf.printme()

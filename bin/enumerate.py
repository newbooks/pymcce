#!/usr/bin/env python
"""
Enumerate the combinations of multiple arrays, one from each.
Sample code to loop over all states in analytical solution.
"""

import itertools

free_residues = [[0,1],[2,3,4],[5,6],[7,8,9]]

old_state = set([x[0] for x in free_residues])
for state in itertools.product(*free_residues):
    print(state)
    new_state = set(state)
    off_conf = list(old_state - new_state)
    on_conf = list(new_state - old_state)
    print("%s -> %s" % (",".join(["%d" % x for x in off_conf]), ",".join(["%d" % x for x in on_conf])))
    old_state = new_state

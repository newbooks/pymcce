#!/usr/bin/env python
import random
import sys
import math
KCAL2KT = 1.688

state_energy = [0.0, 0.0, 0.0]

b = -KCAL2KT


def new_state(state):
    while True:
        try_state = random.randrange(len(state_energy))
        #if try_state != state:
        #    break
        break
    return try_state


def accept(old_state, state):
    dE = state_energy[state] - state_energy[old_state]
    accp = False
    if dE < 0.0:
        accp = True
    elif random.random() < math.exp(b*dE):
        accp = True
    else:
        accp = False
    return accp


if __name__ == "__main__":

    state = random.randrange(len(state_energy))
    occ = [0 for i in range(len(state_energy))]

    if len(sys.argv) < 2:
        print("Specify a sampling number.")
        sys.exit()

    n = int(sys.argv[1])
    for i in range(n):
        old_state = state
        state = new_state(state)
        if not accept(old_state, state):
            state = old_state
        occ[state] += 1

    print(occ)

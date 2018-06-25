#!/usr/bin/env python
import sys


class Env:
    def __init__(self):
        # Hard define values
        self.runprm = "run.prm"
        self.var = {}
        self.param = {}
        return

    def set(self, key, value):
        # Known non-string value are converted, otherwise string presumed
        float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD"
                        "CLASH_DISTANCE"]
        int_values = ["TITR_STEPS"]

        if key in float_values:
            self.var[key] = float(value)
        elif key in int_values:
            self.var[key] = int(value)
        else:
            self.var[key] = value
        return

    def load_runprm(self):
        lines = open(self.runprm).readlines()
        # Sample line: "t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
        for line in lines:
            line = line.strip()
            line = line.split("#")[0]   # This cuts off everything after #
            left_p = line.find("(")
            right_p = line.find(")")
            if left_p > 0 and right_p > left_p+1:
                key = line[left_p + 1:right_p]
                fields = line[:left_p].split()
                if len(fields) < 2:
                    continue
                else:
                    value = fields[0]
            else:
                continue

            self.set(key, value)

        return

    def print_runprm(self):
        for key in self.var.keys():
            print("%-25s:%s" % (key, str(self.var[key])))
        return

    def load_tpl(self, file):
        """Load a tpl file."""
        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            

        return

    def read_extra(self):
        """Read extra.tpl."""
        self.load_tpl(self.var["EXTRA"])
        return

class Protein:
    def __init__(self):
        # Hard coded variables
        self.head3list = "head3.lst"
        return

    def monte_read_energies(self):
        """Read head3.lst, energies, and scale by parameters in extra.tpl."""
        open(self.head3list).readlines()

        return


class Residue:
    def __init__(self):
        return


if __name__ == "__main__":
    env = Env()
    env.load_runprm()
    prot = Protein()
    prot.monte_read_energies()


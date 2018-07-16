#!/usr/bin/env python
import sys, time

int_values = ["TITR_STEPS"]


class Env:
    def __init__(self):
        # Hard define values
        self.runprm = "run.prm"
        self.var = {}
        self.param = {}
        return

    def set(self, key, value):
        # Known non-string value are converted, otherwise string presumed
        float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD", "CLASH_DISTANCE"]
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
                if len(fields) < 1:
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
        float_values = ["EXTRA", "SCALING"]
        int_values = []

        lines = open(file).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys = tuple(keys)

            value_string = fields[1].strip()
            if keys[0] in float_values:
                self.param[keys] = float(value_string)
            elif keys[0] in int_values:
                self.param[keys] = int(value_string)
            else:
                self.param[keys] = value_string

        return

    def read_extra(self):
        """Read extra.tpl."""
        self.load_tpl(self.var["EXTRA"])
        default_values_keys = [("SCALING", "VDW0"),
                               ("SCALING", "VDW1"),
                               ("SCALING", "VDW"),
                               ("SCALING", "TORS"),
                               ("SCALING", "ELE"),
                               ("SCALING", "DSOLV")]
        for element in default_values_keys:
            if element not in self.param:
                self.param[element] = 1.0

        return

    def print_scaling(self):
        """Print scaling factors."""
        #print self.param
        print("   Scaling factors:")
        print("   VDW0  = %.3f" % self.param[("SCALING", "VDW0")])
        print("   VDW1  = %.3f" % self.param[("SCALING", "VDW1")])
        print("   VDW   = %.3f" % self.param[("SCALING", "VDW")])
        print("   TORS  = %.3f" % self.param[("SCALING", "TORS")])
        print("   ELE   = %.3f" % self.param[("SCALING", "ELE")])
        print("   DSOLV = %.3f" % self.param[("SCALING", "DSOLV")])
        print("   Done\n")
        return

class Head3List:
    def __init__(self, fname):
        lines = open(fname).readlines()
        lines.pop(0)
        self.conf_crg = []
        self.conf_name = []
        for line in lines:
            fields = line.split()
            self.conf_name.append(fields[1])
            self.conf_crg.append(float(fields[4]))
        return

class Protein:
    def __init__(self):
        # Hard coded variables
        self.head3list_fname = "head3.lst"
        return

    def monte_read_headlist(self):
        """Read head3.lst."""
        print("   Load conformer list from file %s..." % self.head3list_fname)
        self.head3list = Head3List(self.head3list_fname)

        return


class Residue:
    def __init__(self):
        return


def monte(env):
    """ Monte Carlo sampling """

    timerA = time.time()

    print("   Load environment")
    env.print_scaling()

    prot = Protein()
    prot.monte_read_headlist()




    return


if __name__ == "__main__":
    env = Env()
    env.load_runprm()
    env.read_extra()
    env.print_runprm()

    print("Doing step 4. Monte Carlo sampling")
    monte(env)
    print("End of step 4.")
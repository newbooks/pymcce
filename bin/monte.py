#!/usr/bin/env python
import sys
import time
import os.path

int_values = ["TITR_STEPS"]


class Env:
    def __init__(self):
        # Hard define values
        self.runprm = "run.prm"
        self.constants = {
            "VERSION": "PyMCCE 0.1",
            "FN_CONFLIST1": "head1.lst",
            "FN_CONFLIST2": "head2.lst",
            "FN_CONFLIST3": "head3.lst",
            "ENERGY_TABLE": "energies"
        }
        self.var = {}
        self.param = {}
        self.load_runprm()
        self.read_extra()
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

    def print_param(self):
        for key in self.param.keys():
            print("%-25s:%s" % (key, str(self.param[key])))
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


class Conformer:
    def __init__(self, fields):
        self.confname = fields[1]
        self.flag = fields[2]
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = float(fields[5])
        self.pk0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9])
        self.vdw1 = float(fields[10])
        self.tors = float(fields[11])
        self.epol = float(fields[12])
        self.dsolv = float(fields[13])
        self.extra = float(fields[14])
        self.history = fields[15]
        self.entropy = 0.0   # -TS, will be calculated at entropy sampling
        return


class Protein:
    def __init__(self):
        # Hard coded variables
        self.head3list = []
        self.pairwise = {}
        return

    def load_energy(self, env):
        self.read_headlist(env)
        self.read_pairwise(env)

    def read_headlist(self, env):
        """Read head3.lst."""
        fname = env.constants["FN_CONFLIST3"]
        print("   Load conformer list from file %s..." % fname)
        lines = open(fname).readlines()
        lines.pop(0)
        for line in lines:
            fields = line.split()
            if len(fields) >= 16:
                conf = Conformer(fields)
                self.head3list.append(conf)

        # validate
        confnames = [x.confname for x in self.head3list]
        for name in confnames:
            if len(name) != 14:
                print("%s is not a conformer name.")
                sys.exit()
            occurrence = confnames.count(name)
            if occurrence > 1:
                print("Conformer %s occurred %d times" % occurrence)
                sys.exit()
        return

    def print_headlist(self):
        c = 0
        for conf in self.head3list:
            print("%05d %s" % (c, conf.confname))
            c += 1
        return

    def read_pairwise(self, env):
        """Read pairwise interactions from opp files in folder."""
        folder = env.constants["ENERGY_TABLE"]
        print("   Load pairwise interactions from opp file in folder %s ..." % folder)
        for i in range(len(self.head3list)):
            conf = self.head3list[i]
            oppfile = "%s/%s.opp" % (folder, conf.confname)
            if os.path.isfile(oppfile):
                lines = open(oppfile)
            else:


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
    prot.load_energy(env)
    #prot.print_headlist()



    timerB = time.time()
    print("Total time on MC: %1d seconds.\n" % (timerB-timerA))
    return


if __name__ == "__main__":
    env = Env()

    print("Doing step 4. Monte Carlo sampling")
    monte(env)
    print("End of step 4.")
#!/usr/bin/env python
"""
Global Optimization by Local Equilibrium Sampling
"""
import sys
import os.path
#import numpy as np

ROOMT = 298.15
PH2KCAL = 1.364
pairwise = {}

float_values = ["(EPSILON_PROT)", "(TITR_PH0)", "(TITR_PHD)", "(TITR_EH0)", "(TITR_EHD)", "EXTRA", "SCALING"]
int_values = ["(TITR_STEPS)"]


class Env:
    def __init__(self):
        # Hard define values
        self.fn_runprm = "run.prm"
        self.fn_conflist3 = "head3.lst"
        self.energy_table = "energies"
        self.var = {}

        # load parameters
        self.load_runprm()
        self.read_extra()
        return

    def load_runprm(self):
        lines = open(self.fn_runprm).readlines()
        # Sample line: "t        step 1: pre-run, pdb-> mcce pdb                    (DO_PREMCCE)"
        for line in lines:
            line = line.strip()
            line = line.split("#")[0]  # This cuts off everything after #
            left_p = line.rfind("(")
            right_p = line.rfind(")")

            if left_p > 0 and right_p > left_p + 1:
                key = line[left_p:right_p + 1]
                fields = line[:left_p].split()
                if len(fields) >= 1:
                    value = fields[0]
                    self.set(key, value)
        return

    def set(self, key, value):
        # Known non-string value are converted, otherwise string presumed
        if key in float_values:
            self.var[key] = float(value)
        elif key in int_values:
            self.var[key] = int(value)
        else:
            self.var[key] = value
        return

    def load_tpl(self, fname):
        """Load a tpl file."""
        lines = open(fname).readlines()
        for line in lines:
            line = line.split("#")[0]
            fields = line.split(":")
            if len(fields) != 2:
                continue

            key_string = fields[0].strip()
            keys = key_string.split(",")
            keys = [x.strip().strip("\"") for x in keys]
            keys = [x for x in keys if x]
            keys_str = ",".join(keys)

            value_str = fields[1].strip()
            if keys[0] in float_values:
                self.var[keys_str] = float(value_str)
            elif keys[0] in int_values:
                self.var[keys_str] = int(value_str)
            else:
                self.var[keys_str] = value_str
        return

    def read_extra(self):
        """Read extra.tpl."""
        self.load_tpl(self.var["(EXTRA)"])
        default_values_keys = ["SCALING,VDW0",
                               "SCALING,VDW1",
                               "SCALING,VDW",
                               "SCALING,TORS",
                               "SCALING,ELE",
                               "SCALING,DSOLV"]
        for element in default_values_keys:
            if element not in self.var:
                self.var[element] = 1.0

        return

    def printme(self):
        for key in self.var.keys():
            print("%-25s:%s" % (key, str(self.var[key])))
        return


class Head3lst:
    """ Conformer structure """

    def __init__(self, fields):
        # from head3.lst
        self.iConf = int(fields[0])
        self.confname = fields[1]
        self.flag = fields[2].lower()
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = float(fields[5])
        self.pk0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9]) * env.var["SCALING,VDW0"]
        self.vdw1 = float(fields[10]) * env.var["SCALING,VDW1"]
        self.tors = float(fields[11]) * env.var["SCALING,TORS"]
        self.epol = float(fields[12]) * env.var["SCALING,ELE"]
        self.dsolv = float(fields[13]) * env.var["SCALING,DSOLV"]
        self.extra = float(fields[14])
        self.history = fields[15]
        return

    def printme(self):
        print("%05d %s %c %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s" % (self.iConf,
                                                                                                   self.confname,
                                                                                                   self.flag,
                                                                                                   self.occ,
                                                                                                   self.crg,
                                                                                                   self.em0,
                                                                                                   self.pk0,
                                                                                                   self.ne,
                                                                                                   self.nh,
                                                                                                   self.vdw0,
                                                                                                   self.vdw1,
                                                                                                   self.tors,
                                                                                                   self.epol,
                                                                                                   self.dsolv,
                                                                                                   self.extra,
                                                                                                   self.history))


def load_head3lst():
    conformers = []
    fname = env.fn_conflist3
    lines = open(fname).readlines()
    lines.pop(0)
    for line in lines:
        fields = line.split()
        if len(fields) >= 16:
            conf = Head3lst(fields)
            conformers.append(conf)

    confnames = [x.confname for x in conformers]
    for name in confnames:
        if len(name) != 14:
            print("%s is not a conformer name.")
            sys.exit()
        occurrence = confnames.count(name)
        if occurrence > 1:
            print("Conformer %s occurred %d times" % occurrence)
            sys.exit()
    return conformers


def get_pw(ic, jc):
    if (ic, jc) in pairwise:
        value = pairwise[(ic, jc)]
    else:
        value = 0.0
    return value


def load_pairwise():
    folder = env.energy_table
    n_conf = len(head3lst)

    confnames = [x.confname for x in head3lst]
    scale_ele = env.var["SCALING,ELE"]
    scale_vdw = env.var["SCALING,VDW"]
    for ic in range(n_conf):
        conf = head3lst[ic]
        oppfile = "%s/%s.opp" % (folder, conf.confname)
        if os.path.isfile(oppfile):
            lines = open(oppfile)
            for line in lines:
                fields = line.split()
                if len(fields) < 6:
                    continue
                confname = fields[1]
                jc = confnames.index(confname)
                if jc < 0:
                    print("      Warning: %s in file %s is not a conformer" % (confname, oppfile))
                    continue
                ele = float(fields[2])
                vdw = float(fields[3])
                pairwise[(ic, jc)] = ele * scale_ele + vdw * scale_vdw

    # average the opposite sides
    for ic in range(n_conf-1):
        for jc in range(ic+1, n_conf):
            #if abs(pw[ic][jc] - pw[jc][ic]) > 0.000001:
            #    print("%s %.3f <-> %s %.3f" % (confnames[ic], pw[ic][jc], confnames[jc], pw[jc][ic]))
            averaged_pw = (get_pw(ic, jc) + get_pw(jc, ic)) * 0.5
            if abs(averaged_pw) > 0.0001:
                pairwise[(ic, jc)] = pairwise[(jc, ic)] = averaged_pw
    return

env = Env()
head3lst = load_head3lst()
load_pairwise()

#if __name__ == "__main__":

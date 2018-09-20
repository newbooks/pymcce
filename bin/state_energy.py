#!/usr/bin/env python
import sys
import time
import os.path
from copy import deepcopy
import random

ROOMT = 298.15
PH2KCAL = 1.364

class Env:
    def __init__(self):
        # Hard define values
        self.runprm = "run.prm"
        self.version = "PyMCCE 0.1"
        self.fn_conflist1 = "head1.lst"
        self.fn_conflist2 = "head2.lst"
        self.fn_conflist3 = "head3.lst"
        self.energy_table = "energies"
        # run.prm parameters key:value
        self.var = {}
        # tpl parameters (key1, key2, key3):value
        self.param = {}

        # load parameters
        self.load_runprm()
        self.read_extra()
        return

    def set(self, key, value):
        # Known non-string value are converted, otherwise string presumed
        float_values = ["EPSILON_PROT", "TITR_PH0", "TITR_PHD", "TITR_EH0", "TITR_EHD", "CLASH_DISTANCE",
                        "BIG_PAIRWISE", "MONTE_T"]
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
            line = line.split("#")[0]  # This cuts off everything after #
            left_p = line.find("(")
            right_p = line.find(")")

            if left_p > 0 and right_p > left_p + 1:
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
        # print self.param
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
        self.iConf = int(fields[0])
        self.confname = fields[1]
        self.flag = fields[2]
        self.occ = float(fields[3])
        self.crg = float(fields[4])
        self.em0 = float(fields[5])
        self.pk0 = float(fields[6])
        self.ne = int(fields[7])
        self.nh = int(fields[8])
        self.vdw0 = float(fields[9]) * env.param[("SCALING", "VDW0")]
        self.vdw1 = float(fields[10]) * env.param[("SCALING", "VDW1")]
        self.tors = float(fields[11]) * env.param[("SCALING", "TORS")]
        self.epol = float(fields[12]) * env.param[("SCALING", "ELE")]
        self.dsolv = float(fields[13]) * env.param[("SCALING", "DSOLV")]
        self.extra = float(fields[14])
        self.history = fields[15]
        self.entropy = 0.0  # -TS, will be calculated at entropy sampling
        return

    def printme(self):
        print("%05d %s %c %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s" % (self.iConf,
                                                                                                   self.confname,
                                                                                                   self.flag, self.occ,
                                                                                                   self.crg,
                                                                                                   self.em0, self.pk0,
                                                                                                   self.ne, self.nh,
                                                                                                   self.vdw0, self.vdw1,
                                                                                                   self.tors, self.epol,
                                                                                                   self.dsolv,
                                                                                                   self.extra,
                                                                                                   self.history))


class Protein:
    def __init__(self):
        # Hard coded variables
        self.head3list = []
        self.confnames = []
        self.pairwise = {}
        self.fixed_conformers = []
        self.free_residues = []
        self.self_energy = []
        self.ph = 0.0
        self.eh = 0.0
        return

    def load_energy(self):
        self.read_headlist()
        self.read_pairwise()

    def read_headlist(self):
        """Read head3.lst."""
        fname = env.fn_conflist3
        print("   Load conformer list from file %s..." % fname)
        lines = open(fname).readlines()
        lines.pop(0)
        for line in lines:
            fields = line.split()
            if len(fields) >= 16:
                conf = Conformer(fields)
                self.head3list.append(conf)

        # validate
        self.confnames = [x.confname for x in self.head3list]
        for name in self.confnames:
            if len(name) != 14:
                print("%s is not a conformer name.")
                sys.exit()
            occurrence = self.confnames.count(name)
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

    def read_pairwise(self):
        """Read pairwise interactions from opp files in folder."""
        folder = env.energy_table
        print("   Load pairwise interactions from opp file in folder %s ..." % folder)
        for i in range(len(self.head3list)):
            conf = self.head3list[i]
            oppfile = "%s/%s.opp" % (folder, conf.confname)
            resid_i = conf.confname[:3] + conf.confname[5:11]
            if os.path.isfile(oppfile):
                lines = open(oppfile)
                for line in lines:
                    fields = line.split()
                    if len(fields) < 6:
                        continue
                    confname = fields[1]
                    j = self.confnames.index(confname)
                    if j < 0:
                        print("      Warning: %s in file %s is not a conformer" % (confname, oppfile))
                        continue

                    resid_j = confname[:3] + confname[5:11]
                    if resid_i != resid_j:  # not within a residue
                        ele = float(fields[2])
                        vdw = float(fields[3])
                        self.pairwise[(i, j)] = ele * env.param[("SCALING", "ELE")] + vdw * env.param[("SCALING", "VDW")]
            else:  # No opp files, assume all interactions to be 0, and no entries in the pairwise{}
                pass
        return

    def update_self_energy(self):
        # This changes with entropy sampling
        for conf in self.head3list:
            monte_temp = env.var["MONTE_T"]
            E_ph = monte_temp / ROOMT * conf.nh * (self.ph - conf.pk0) * PH2KCAL
            E_eh = monte_temp / ROOMT * conf.ne * (self.eh - conf.em0) * PH2KCAL/58.0
            Eself = conf.vdw0 + conf.vdw1 + conf.epol + conf.tors + conf.dsolv + conf.extra + E_ph + E_eh + conf.entropy
            self.self_energy.append(Eself)
        return

    def group_to_residues(self):
        # Residues will be put into an irregular 2D array, the rows are residue index and the elements in row are
        # conformer indices.
        # Residues will be further divided into fixed residues and free residues upon examination
        residue_ids = []
        for confname in self.confnames:
            resid = confname[:3] + confname[5:11]
            if resid not in residue_ids:
                residue_ids.append(resid)
        self.residues = [[] for i in range(len(residue_ids))]
        for i in range(len(self.confnames)):
            confname = self.confnames[i]
            resid = confname[:3] + confname[5:11]
            index = residue_ids.index(resid)
            self.residues[index].append(i)

        # Verify head3list flag and occ; Find free and fixed residues
        # if total occ of "t" flagged conformers is 1:
        # the rest conformers will be set to "t 0.00", and this residue is "fixed"
        # else if total occ of "t" flagged conformers is 0:
        #    if only one conformer is left:
        #        the lone conformer is set to "t 1.00", and this conformer and this residue will be "fixed"
        #    else:
        #        this residue is "free" and occ of conformers is 0.
        # otherwise:
        #    partial fixed occupancy not allowed

        print("   Grouping and verifying conformers ...")
        for res in self.residues:
            socc = 0.0
            n_freeconf = len(res)
            for i in res:
                if self.head3list[i].flag.upper() == "T":
                    socc += self.head3list[i].occ
                    n_freeconf -= 1
            if abs(socc - 1.0) < 0.001:  # total occ of fixed conformers are 1.0
                for i in res:
                    self.fixed_conformers.append(i)
                    if self.head3list[i].flag.upper() != "T":
                        print("      %s %c %4.2f -> %s t  0.00" % (self.head3list[i].confname, self.head3list[i].flag, self.head3list[i].occ, self.head3list[i].confname))
                        self.head3list[i].occ = 0.0
                        self.head3list[i].flag = "t"
            elif abs(socc) < 0.001:  # total occ is 0
                if n_freeconf == 1:
                    for i in res:
                        if self.head3list[i].flag.upper() != "T":
                            print("      %s %c %4.2f -> %s t  1.00" % (self.head3list[i].confname, self.head3list[i].flag, self.head3list[i].occ, self.head3list[i].confname))
                            self.head3list[i].flag = "t"
                            self.head3list[i].occ = 1.0
                            self.fixed_conformers.append(i)
                            break  # because only one "f"
                else:
                    free_conformers = []
                    for i in res:
                        if self.head3list[i].flag.upper() == "T":
                            self.fixed_conformers.append(i)
                        else:
                            free_conformers.append(i)
                    self.free_residues.append(free_conformers)
            else:  # total occ is neither 0 or 1
                print("   Error: Total residue occupancy is %.2f, 0.00 or 1.00 expected." % socc)
                for i in res:
                    self.head3list[i].printme()
                print("   Exiting ...")
                sys.exit()

        # Make a running copy and update biglist
        for conf in self.head3list:
            conf.fixed_occ = conf.occ
        self.fixed_conformers_running = deepcopy(self.fixed_conformers)
        self.free_residues_running = deepcopy(self.free_residues)
        self.biglist = self.make_biglist(self.free_residues_running)
        return

    def make_biglist(self, free_residues):
        biglist = [[] for i in range(len(free_residues))]
        for ir in range(len(free_residues)):
            for jr in range(ir+1, len(free_residues)):
                next_jr = False
                for ic in free_residues[ir]:
                    if next_jr:
                        break
                    for jc in free_residues[jr]:
                        if next_jr:
                            break
                        if (ic, jc) in self.pairwise:
                            pw = self.pairwise[(ic, jc)]
                        else:
                            pw = 0.0
                        if abs(pw) > env.var["BIG_PAIRWISE"]:
                            biglist[ir].append(jr)
                            biglist[jr].append(ir)
                            next_jr = True
        return biglist


def randomize_state(free_res):
    state = []
    for res in free_res:
        if len(res) < 2:
            print("   Error: Randomize a residue with less than 2 conformers")
            sys.exit()

        state.append(res[random.randrange(len(res))])
    return state


def get_E(complete_state, prot):
    # This function requires a complete state and supports partial occupancy
    E_state = 0.0
    for ic in complete_state:
        if ic in prot.fixed_conformers:
            E_state += prot.self_energy[ic] * prot.head3list[ic].fixed_occ
        else:
            E_state += prot.self_energy[ic]

    # To support partial occ, multiply fixed_occ for fixed conformers
    for i in range(len(complete_state) - 1):
        for j in range(i+1, len(complete_state)):
            ic = complete_state[i]
            jc = complete_state[j]
            pair = (ic, jc)
            #print pair, prot.pairwise[pair]
            if pair in prot.pairwise:
                pw_energy = prot.pairwise[pair]
                if ic in prot.fixed_conformers_running:
                    pw_energy *= prot.head3list[ic].fixed_occ
                if jc in prot.fixed_conformers_running:
                    pw_energy *= prot.head3list[jc].fixed_occ
                E_state += pw_energy

    return E_state


def monte():
    """ Monte Carlo sampling """

    timerA = time.time()

    print("   Load environment")
    env.print_scaling()

    prot = Protein()
    prot.load_energy()
    prot.group_to_residues()

    # Gnerate a microstate in size of prot.free_residues_running[]
    state = randomize_state(prot.free_residues_running)
    complete_state = prot.fixed_conformers_running + state
    complete_state.sort()

    #print complete_state

    prot.ph = 7.0
    prot.eh = 0.0
    prot.update_self_energy()

    lines = open("states").readlines()
    for line in lines:
        state = [int(x.strip(",")) for x in line.split(",")]
        E_state = get_E(state, prot)
        print("[%s] E_state=%.4f" % (" ".join(str(x) for x in state), E_state))


    #E_state = get_E(complete_state, prot)

    #print(complete_state)
    #print("E_state = %.2f" % E_state)

    timerB = time.time()
    print("Total time on MC: %1d seconds.\n" % (timerB - timerA))
    return


if __name__ == "__main__":
    env = Env()

    print("Doing step 4. Monte Carlo sampling")
    monte()
    print("End of step 4.")

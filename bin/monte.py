#!/usr/bin/env python
import sys
import time
import os.path
from copy import deepcopy
import random
import math
import numpy as np

ROOMT = 298.15
PH2KCAL = 1.364
KCAL2KT = 1.688
KJ2KCAL = 0.239


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
                        "BIG_PAIRWISE", "MONTE_T", "MONTE_REDUCE"]
        int_values = ["TITR_STEPS", "MONTE_RUNS", "MONTE_TRACE", "MONTE_NITER", "MONTE_NEQ",
                      "MONTE_NSTART", "MONTE_FLIPS"]
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
            left_p = line.rfind("(")
            right_p = line.rfind(")")

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
        self.vdw0 = float(fields[9]) * env.param[("SCALING", "VDW0")]
        self.vdw1 = float(fields[10]) * env.param[("SCALING", "VDW1")]
        self.tors = float(fields[11]) * env.param[("SCALING", "TORS")]
        self.epol = float(fields[12]) * env.param[("SCALING", "ELE")]
        self.dsolv = float(fields[13]) * env.param[("SCALING", "DSOLV")]
        self.extra = float(fields[14])
        self.history = fields[15]
        # from MC entropy sampling
        self.entropy = 0.0  # -TS, will be calculated at entropy sampling
        # needed by MC process
        self.E_self = 0.0  # self energy in head3.lst
        self.E_self_mfe = 0.0  # self energy including pairwise contribution from fixed residues
        self.counter = 0  # MC counters
        self.mc_occ = 0.0  # MC calculated occ
        self.acc_occ = []  # a list of past mc_occ until rest, history is needed for test convergence
        self.occ_at_points = []  # occ at titration points
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
        self.ph = 0.0
        self.eh = 0.0
        return

    def reset_flag(self):
        for conf in self.head3list:
            if conf.flag == "t":  # fixed occ, not titratable
                conf.on = False
            else:
                conf.on = True
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

    def print_entropy(self):
        c = 0
        for conf in self.head3list:
            mc_log.write("   %05d %s mc_occ = %4.2f TS = %5.2f\n" % (c, conf.confname, conf.mc_occ, conf.entropy))
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
                        pw = ele * env.param[("SCALING", "ELE")] + vdw * env.param[("SCALING", "VDW")]
                        if (j, i) in self.pairwise:  # average if other direction exists
                            # print("%.2f %.2f" % (pw, self.pairwise[(j, i)]))
                            pw = 0.5 * (pw + self.pairwise[(j, i)])
                        self.pairwise[(i, j)] = self.pairwise[
                            (j, i)] = pw  # use the same value for both directions by default

            else:  # No opp files, assume all interactions to be 0, and no entries in the pairwise{}
                pass
        return


    def verify_group_conformers(self):
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
        # if only one conformer is left:
        #        the lone conformer is set to "t 1.00", and this conformer and this residue will be "fixed"
        #    else:
        #        this residue is "free" and occ of conformers is 0.
        # otherwise:
        #    partial fixed occupancy not allowed

        print("   Verifying conformers ...")
        # Verify flags
        # Group conformers
        for res in self.residues:
            socc = 0.0
            n_freeconf = len(res)
            for i in res:
                if not self.head3list[i].on:
                    socc += self.head3list[i].occ
                    n_freeconf -= 1
                elif abs(self.head3list[i].occ) > 0.001:  # free residue has non-0 occ
                    print("      %s %c %4.2f -> %s f  0.00 (free conformer initial occ = 0)" % (
                        self.head3list[i].confname,
                        self.head3list[i].flag,
                        self.head3list[i].occ, self.head3list[i].confname))
                    self.head3list[i].occ = 0.0
            if abs(socc - 1.0) < 0.001:  # total occ of fixed conformers are 1.0
                for i in res:
                    self.fixed_conformers.append(i)
                    if self.head3list[i].on:
                        print("      %s %c %4.2f -> %s t  0.00 (fixed conformers already have occ 1.0)" % (
                            self.head3list[
                                i].confname,
                            self.head3list[i].flag, self.head3list[i].occ, self.head3list[i].confname))
                        self.head3list[i].occ = 0.0
                        self.head3list[i].on = False
                        self.head3list[i].flag = "t"
            elif abs(socc) < 0.001:  # total occ is 0
                if n_freeconf == 1:
                    for i in res:
                        if self.head3list[i].on:
                            print("      %s %c %4.2f -> %s t  1.00 (single conformer of the residue)" % (self.head3list[
                                                                                                             i].confname,
                                                                                                         self.head3list[
                                                                                                             i].flag,
                                                                                                         self.head3list[
                                                                                                             i].occ,
                                                                                                         self.head3list[
                                                                                                             i].confname))
                            self.head3list[i].on = False
                            self.head3list[i].occ = 1.0
                            self.head3list[i].flag = "t"
                            self.fixed_conformers.append(i)
                            break  # because only one "f"
                else:
                    free_conformers = []
                    for i in res:
                        if not self.head3list[i].on:
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

#        print(self.fixed_conformers)
#        for x in self.free_residues:
#            print x

        return

    def group_conformers(self):
        residue_ids = []
        self.fixed_conformers = []
        self.free_residues = []
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

        for res in self.residues:
            socc = 0.0
            n_freeconf = len(res)
            for i in res:
                if not self.head3list[i].on:
                    socc += self.head3list[i].occ
                    n_freeconf -= 1
                elif abs(self.head3list[i].occ) > 0.001:  # free residue has non-0 occ
                    self.head3list[i].occ = 0.0
            if abs(socc - 1.0) < 0.001:  # total occ of fixed conformers are 1.0
                for i in res:
                    self.fixed_conformers.append(i)
                    if self.head3list[i].on:
                        self.head3list[i].occ = 0.0
                        self.head3list[i].on = False
            elif abs(socc) < 0.001:  # total occ is 0
                if n_freeconf == 1:
                    for i in res:
                        if self.head3list[i].on:
                            self.head3list[i].on = False
                            self.head3list[i].occ = 1.0
                            self.fixed_conformers.append(i)
                            break  # because only one "f"
                else:
                    free_conformers = []
                    for i in res:
                        if not self.head3list[i].on:
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

#        print(self.fixed_conformers)
#        for x in self.free_residues:
#            print x

        return


def get_entropy(prot):
    # calculate entropy from mc_occ and return the biggest difference from previous entropy
    max_dTS = 0.0

    for res in prot.residues:
        typeids = {}
        for ic in res:
            # compose typeid from conformrs within a residue
            # this should include ne, nh, and DM status (the 4th character of conformer name)
            typeid = (prot.head3list[ic].ne, prot.head3list[ic].nh, prot.head3list[ic].confname[3])
            if typeid in typeids:
                typeids[typeid].append(ic)
            else:
                typeids[typeid] = [ic]

        for typeid in typeids:
            TS = 0.0
            sum = 0.0
            confs = typeids[typeid]

            for ic in confs:
                if prot.head3list[ic].on:
                    sum += prot.head3list[ic].mc_occ
                else:
                    sum += prot.head3list[ic].occ
                    prot.head3list[ic].mc_occ = prot.head3list[ic].occ

            if sum < 0.0001:
                for ic in confs:
                    if max_dTS < abs(prot.head3list[ic].entropy):
                        max_dTS = abs(prot.head3list[ic].entropy)
                    prot.head3list[ic].entropy = 0.0
            else:
                for ic in confs:
                    p = prot.head3list[ic].mc_occ/sum
                    if p < 1.0E-6:
                        TS = 0.0
                    else:
                        TS -= p * math.log(p) / 1.688

                for ic in confs:
                    if max_dTS < abs(prot.head3list[ic].entropy - TS):
                        max_dTS = abs(prot.head3list[ic].entropy - TS)
                    prot.head3list[ic].entropy = TS

    return max_dTS




class MicroState:
    def __init__(self, prot):
        self.E_state = 0.0
        self.fixed_conformers = prot.fixed_conformers
        self.free_residues = deepcopy(prot.free_residues)
        self.biglist = self.make_biglist(prot)
        self.state = self.randomize_state()
        self.complete_state = []
        return

    def update_self_energy(self, prot):
        # This changes with entropy sampling
        for ic in range(len(prot.head3list)):
            conf = prot.head3list[ic]
            monte_temp = env.var["MONTE_T"]
            E_ph = monte_temp / ROOMT * conf.nh * (prot.ph - conf.pk0) * PH2KCAL
            E_eh = monte_temp / ROOMT * conf.ne * (prot.eh - conf.em0) * PH2KCAL / 58.0
            Eself = conf.vdw0 + conf.vdw1 + conf.epol + conf.tors + conf.dsolv + conf.extra + E_ph + E_eh + \
                    conf.entropy # entropy term is an artificial "correction", thus add
            prot.head3list[ic].E_self = Eself

            # mfe from fixed conformer
            mfe = 0.0
            for jc in self.fixed_conformers:
                if (ic, jc) in prot.pairwise:
                    mfe += prot.pairwise[(ic, jc)] * prot.head3list[jc].occ
            prot.head3list[ic].E_self_mfe = prot.head3list[ic].E_self + mfe

        return

    def randomize_state(self):
        state = []
        for res in self.free_residues:
            if len(res) < 2:
                print("   Error: Randomize a residue with less than 2 conformers")
                sys.exit()
            state.append(res[random.randrange(len(res))])
        return state

    def make_biglist(self, prot):
        biglist = [[] for i in range(len(self.free_residues))]
        for ir in range(len(self.free_residues)):
            for jr in range(ir + 1, len(self.free_residues)):
                next_jr = False
                for ic in self.free_residues[ir]:
                    if next_jr:
                        break
                    for jc in self.free_residues[jr]:
                        if next_jr:
                            break
                        if (ic, jc) in prot.pairwise:
                            pw = prot.pairwise[(ic, jc)]
                        else:
                            pw = 0.0
                        if abs(pw) > env.var["BIG_PAIRWISE"]:
                            biglist[ir].append(jr)
                            biglist[jr].append(ir)
                            next_jr = True
        return biglist


    def get_E(self, prot):
        self.complete_state = self.fixed_conformers + self.state
        # This function requires a complete state and supports partial occupancy
        E_state = 0.0
        for ic in self.complete_state:
            if ic in self.fixed_conformers:
                E_state += prot.head3list[ic].E_self * prot.head3list[ic].occ
            else:
                E_state += prot.head3list[ic].E_self

        # To support partial occ, multiply occ for fixed conformers
        for i in range(len(self.complete_state) - 1):
            for j in range(i + 1, len(self.complete_state)):
                ic = self.complete_state[i]
                jc = self.complete_state[j]
                pair = (ic, jc)
                # print pair, prot.pairwise[pair]
                if pair in prot.pairwise:
                    pw_energy = prot.pairwise[pair]
                    if ic in self.fixed_conformers:
                        pw_energy *= prot.head3list[ic].occ
                    if jc in self.fixed_conformers:
                        pw_energy *= prot.head3list[jc].occ
                    E_state += pw_energy

        return E_state


    def verify_free(self):
        for x in self.free_residues:
            if len(x) < 2:
                print("   ERROR: Free residues can not have less than 2 flippable conformers.")
                sys.exit()
        return


def mc_complete_sampling(prot):
    global mc_log, ms_dat

    prot.reset_flag()
    prot.group_conformers()
    state = MicroState(prot)
    state.update_self_energy(prot)
    state.biglist = state.make_biglist(prot)

    # Annealing`
    mc_log.write("Doing annealing:\n")
    # At this time, fixed_conformers are logically determined
    N = env.var["MONTE_NSTART"] * (len(prot.head3list) - len(prot.fixed_conformers)) / 10
    for i in range(5):
        mc_run(prot, state, N, T=env.var["MONTE_T"] + (9 - i) * 100.0)
    # Run extended time at the lowest temperature
    mc_run(prot, state, N + 1000, T=env.var["MONTE_T"])

    # Reduction
    mc_log.write("Doing reduction:\n")
    N = env.var["MONTE_NEQ"] * (len(prot.head3list) - len(prot.fixed_conformers))
    mc_run(prot, state, N, T=env.var["MONTE_T"], record=True)
    for res_confs in state.free_residues:
        for ic in res_confs:
            if prot.head3list[ic].mc_occ < env.var["MONTE_REDUCE"]:
                prot.head3list[ic].on = False
                prot.head3list[ic].occ = 0.0

    prot.group_conformers()
    state = MicroState(prot)
    state.update_self_energy(prot)
    state.biglist = state.make_biglist(prot)

    # Entropy sampling
    # Use the last reduction run to estimate entropy
    get_entropy(prot)
    state.update_self_energy(prot)

    # Throw away runs because a new random state was created
    N = env.var["MONTE_NEQ"] * (len(prot.head3list) - len(prot.fixed_conformers))
    mc_run(prot, state, N, T=env.var["MONTE_T"])

    early_break = False
    for i in range(6):
        mc_log.write("Doing entropy sampling...\n")
        mc_run(prot, state, N, T=env.var["MONTE_T"], record=True)
        max_dTS = get_entropy(prot)
        state.update_self_energy(prot)

        # prot.print_entropy()
        mc_log.write("   Maximum entropy difference = %.2f\n" % max_dTS)
        if max_dTS <= 0.5:
            early_break = True
            break

    if early_break:
        mc_log.write("   End entropy sampling\n")
    else:
        mc_log.write("   Entropy correction didn't converge within 6 time.\n")

    prot.print_entropy()

    # Sampling
    mc_log.write("Doing sampling:\n")
    N = env.var["MONTE_NITER"] * (len(prot.head3list) - len(prot.fixed_conformers))
    mc_run(prot, state, N, record=True)

    for conf in prot.head3list:
        conf.acc_occ.append(conf.mc_occ)


def mc_run(prot, state, N, T=ROOMT, record=False):
    """Monte Carlo sampling core function. It starts from a state, sample N times."""

    global mc_log
    global ms_dat

    state.verify_free()

    if "MONTE_WRITESTATE" in env.var and env.var["MONTE_WRITESTATE"].upper() == "T":
        write_ms = True
    else:
        write_ms = False

    b = -KCAL2KT / (T / ROOMT)
    n_free = len(state.free_residues)
    nflips = env.var["MONTE_FLIPS"]

    E_minimum = E_state = state.get_E(prot)

    # Trace cycles
    if env.var["MONTE_TRACE"] > 0:
        cycles = int((N - 1) / env.var["MONTE_TRACE"]) + 1  # minimum 1
        n_total = cycles * env.var["MONTE_TRACE"]
        n_cycle = env.var["MONTE_TRACE"]
    else:
        cycles = 1
        n_total = n_cycle = N

    # clear counters
    for conf in prot.head3list:
        conf.counter = 0
    H_average = 0.0


    # clear counters
    for res_confs in state.free_residues:
        for ic in res_confs:
            prot.head3list[ic].counter = 0

    for i in range(cycles):
        mc_log.write("   Step %10d, E_minimum = %10.2f, E_running = %10.2f\n" % (i * n_cycle, E_minimum, E_state))
        mc_log.flush()

        iters = n_cycle
        while iters:
            # save the state
            old_state = [x for x in state.state]

            # 1st flip
            ires = random.randrange(n_free)
            while True:
                new_conf = random.choice(state.free_residues[ires])
                if new_conf != state.state[ires]:
                    break

            old_conf = state.state[ires]
            state.state[ires] = new_conf
            dE = prot.head3list[new_conf].E_self_mfe - prot.head3list[old_conf].E_self_mfe
            for j in range(n_free):
                if (new_conf, state.state[j]) in prot.pairwise:
                    dE += prot.pairwise[(new_conf, state.state[j])]
                if (old_conf, state.state[j]) in prot.pairwise:
                    dE -= prot.pairwise[(old_conf, state.state[j])]

            if random.choice([True, False]):
                if state.biglist[ires]:
                    for k in range(nflips):
                        iflip = random.choice(state.biglist[ires])  # which residue to flip
                        old_conf = state.state[iflip]
                        new_conf = random.choice(state.free_residues[iflip])  # conformer flip to
                        state.state[iflip] = new_conf
                        # print("iflip=%d from biglist=%s ")

                        dE += prot.head3list[new_conf].E_self_mfe - prot.head3list[old_conf].E_self_mfe
                        for j in range(n_free):
                            if (new_conf, state.state[j]) in prot.pairwise:
                                dE += prot.pairwise[(new_conf, state.state[j])]
                            if (old_conf, state.state[j]) in prot.pairwise:
                                dE -= prot.pairwise[(old_conf, state.state[j])]

            if dE < 0.0:
                flip = True
            elif random.random() < math.exp(b * dE):
                flip = True
            else:
                flip = False

            if flip:  # update energy
                E_state += dE
                if E_minimum > E_state:
                    E_minimum = E_state
            else:  # go back to old state
                state.state = [x for x in old_state]

            if record:  # record this to conformer stat
                for ic in state.state:
                    prot.head3list[ic].counter += 1

            H_average += E_state
            iters -= 1

    # find mc_occ
    if record:
        for res_confs in state.free_residues:
            for ic in res_confs:
                prot.head3list[ic].mc_occ = float(prot.head3list[ic].counter) / n_total


    state.complete_state = state.state + state.fixed_conformers
    state.complete_state.sort()
    state_str = ",".join(["%d" % x for x in state.complete_state])
    mc_log.write("%s\n" % state_str)
    mc_log.write("   Exit %10d, E_minimum = %10.2f, E_running = %10.2f\n" % (n_total, E_minimum, E_state))
    mc_log.write("   The average running energy, corresponding to H, is %8.3f\n" % (H_average / n_total))
    mc_log.flush()
    return


def monte():
    """ Monte Carlo sampling """
    global mc_log

    timerA = time.time()

    print("   Load environment")
    env.print_scaling()

    prot = Protein()
    prot.load_energy()
    prot.reset_flag()
    prot.verify_group_conformers()
    prot.monte_t = env.var["MONTE_T"]

    titration_type = env.var["TITR_TYPE"].upper()
    init_ph = env.var["TITR_PH0"]
    step_ph = env.var["TITR_PHD"]
    init_eh = env.var["TITR_EH0"]
    step_eh = env.var["TITR_EHD"]
    steps = env.var["TITR_STEPS"]
    mc_log = open("mc.out", "w")

    timerB = time.time()
    print("   Done setting up MC in %1d seconds.\n" % (timerB - timerA))

    points = []
    for i in range(steps):
        # Set up pH and eh environment
        if titration_type == "PH":
            prot.ph = init_ph + i * step_ph
            prot.eh = init_eh
            points.append(prot.ph)
        elif titration_type == "EH":
            prot.ph = init_ph
            prot.eh = init_eh + i * step_eh
            points.append(prot.eh)
        else:
            print(
            "   Error: Titration type is %s. It has to be ph or eh in line (TITR_TYPE) in run.prm" % titration_type)
            sys.exit()

        print("      Titration at ph = %5.2f and eh = %.0f mv." % (prot.ph, prot.eh))
        mc_log.write("Titration at ph = %5.2f and eh = %.0f mv.\n" % (prot.ph, prot.eh))
        # Gnerate a microstate in size of prot.free_residues_running[]


        # Reset mc_occ accumulated list
        for ic in range(len(prot.head3list)):
            prot.head3list[ic].acc_occ = []

        for j_runs in range(env.var["MONTE_RUNS"]):
            # independent runs
            mc_complete_sampling(prot)

        mc_log.write("Standard deviation in independent runs:\n")
        max_std = 0.0
        for conf in prot.head3list:
            occs = ["%4.2f" % x for x in conf.acc_occ]
            a = np.array(conf.acc_occ)
            s = np.std(a)
            m = np.mean(a)
            conf.occ_at_points.append(float(m))
            if s > max_std:
                max_std = s
            mc_log.write("   %s %s stdev = %6.3f\n" % (conf.confname, " ".join(occs), float(s)))

        if max_std > 0.1:
            print("      Warning: found a large standard deveiation %.3f. Check mc.log for details" % max_std)



    # print fort.38
    lines = []
    if titration_type.upper == "PH":
        points_str = " ".join(["%5.3f" % x for x in points])
    else:
        points_str = " ".join(["%5.f" % x for x in points])
    lines.append("        %s     %s\n" % (titration_type, points_str))

    for conf in prot.head3list:
        occ_str = " ".join(["%5.3f" % x for x in conf.occ_at_points])
        line = "%s %s\n" % (conf.confname, occ_str)
        lines.append(line)

    open("fort.38", "w").writelines(lines)


    mc_log.close()
    timerB = time.time()
    print("Total time on MC: %1d seconds.\n" % (timerB - timerA))
    return


if __name__ == "__main__":
    env = Env()

    print("Doing step 4. Monte Carlo sampling")
    monte()
    print("End of step 4.")

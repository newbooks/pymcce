#!/usr/bin/env python
"""
Global Optimization by Local Equilibrium Sampling
"""
import sys
import os.path
import random
import numpy as np

ROOMT = 298.15
PH2KCAL = 1.364
CLUSTER_PWCUTOFF = 0.5  # include into a cluster if conf-conf pw is bigger than this value
CLUSTER_EXTENDED_LEVEL = 3  # extend nodes at this level of neighbor search
EQ_MAXCYCLES = 100   # maximum cycles to equilibrate the clusters
EQ_CONVERGE = 0.01   # occupancy convergence (standard deviation of free confs) criterion


residue_report = "nodes.info"
neighbor_report = "neighbor.info"
cluster_report = "cluster.info"
converge_progress = "converge.progress"
energy_progress = "energy.progress"
mc_progress = "mc.log"

float_values = ["(EPSILON_PROT)", "(TITR_PH0)", "(TITR_PHD)", "(TITR_EH0)", "(TITR_EHD)", "EXTRA", "SCALING",
                "(MONTE_T)"]
int_values = ["(TITR_STEPS)", "(NSTATE_MAX)"]


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


class Residue:
    def __init__(self, resid):
        self.resid = resid
        self.neighbors = []
        self.conformers = []
        self.fixed_conformers = []
        self.free_conformers = []
        return

    def verify_flags(self):
        socc = 0.0
        n_freeconf = len(self.conformers)
        for conf in self.conformers:
            if not conf.on:
                socc += conf.occ
                n_freeconf -= 1
            elif abs(conf.occ) > 0.001:  # free residue has non-0 occ
                print("      %s %c %4.2f -> %s f  0.00 (free conformer initial occ = 0)" % (
                    head3lst[conf.i].confname,
                    conf.flag,
                    conf.occ, head3lst[conf.i].confname))
                conf.occ = 0.0
        if abs(socc - 1.0) < 0.001:  # total occ of fixed conformers are 1.0, all fixed
            for conf in self.conformers:
                if conf.on:
                    print("      %s %c %4.2f -> %s t  0.00 (fixed conformers already have occ 1.0)" % (
                        head3lst[conf.i].confname,
                        conf.flag, conf.occ, head3lst[conf.i].confname))
                    conf.occ = 0.0
                    conf.on = False
                    conf.flag = "t"
                self.fixed_conformers.append(conf)
        elif abs(socc) < 0.001:  # total occ is 0
            if n_freeconf == 1:  # The only free conformer will be rest to t 1.0
                for conf in self.conformers:
                    if conf.on:
                        print("      %s %c %4.2f -> %s t  1.00 (single free conformer of the residue)" % (
                            head3lst[conf.i].confname,
                                                                                                          conf.flag,
                                                                                                          conf.occ,
                                                                                                          head3lst[conf.i].confname))
                        conf.on = False
                        conf.occ = 1.0
                        conf.flag = "t"
                        self.fixed_conformers.append(conf)
                    else:
                        self.fixed_conformers.append(conf)  # The rest are fixed at 0
            else:
                for conf in self.conformers:
                    if not conf.on:
                        self.fixed_conformers.append(conf)
                    else:
                        self.free_conformers.append(conf)
        else:  # total occ is neither 0 or 1
            print("   Error: Total residue occupancy is %.2f, 0.00 or 1.00 expected." % socc)
            for conf in self.conformers:
                conf.printme()
            print("   Exiting ...")
            sys.exit()
        return


class Conformer:
    def __init__(self, ic):
        self.i = ic    # make a link to the head3lst
        self.flag = ""
        self.on = True  # True means to be sampled, False means fixed at occ
        self.occ = 0.0
        self.occ_old = 0.0
        self.counter = 0
        self.E_self = head3lst[ic].vdw0 + head3lst[ic].vdw1 + head3lst[ic].epol + head3lst[ic].tors + head3lst[ic].dsolv + head3lst[ic].extra
        self.E_pheh = 0.0  # ph and eh effect
        self.E_mfe = 0.0   # when in an active cluster, this is the pairwise from outside conformers.
        self.entropy = 0.0 # entropy correction, used by metropolis criterion but not by system energy.
        self.entropy_old = 0.0
        self.E_total = 0.0
        self.type = ""
        self.load(ic)
        return

    def load(self, ic):
        self.flag = head3lst[ic].flag
        if head3lst[ic].flag.upper() == "T":
            self.on = False
            self.occ = head3lst[ic].occ
        else:
            self.on = True
            self.occ = 0.0
        return

    def printme(self):
        if self.on:
            flag = "f"
        else:
            flag = "t"
        line = "%05d %s %c %4.2f %6.3f %5d %5.2f %2d %2d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s" % (head3lst[self.i].iConf,
               head3lst[self.i].confname,
               flag,
               self.occ,
               head3lst[self.i].crg,
               head3lst[self.i].em0,
               head3lst[self.i].pk0,
               head3lst[self.i].ne,
               head3lst[self.i].nh,
               head3lst[self.i].vdw0,
               head3lst[self.i].vdw1,
               head3lst[self.i].tors,
               head3lst[self.i].epol,
               head3lst[self.i].dsolv,
               head3lst[self.i].extra,
               head3lst[self.i].history)
        return line


class Node:
    """
    Node is a simplified structure for free residues and contains only free conformers. It is the core for MC sampling.
    """
    def __init__(self, res):
        self.resid = res.resid
        self.free_conformers = res.free_conformers
        self.neighbors = []
        return

    def find_neighbors(self):
        for n in nodes:
            if self.resid != n.resid:
                maxpw = self.find_maxpw(self, n)
                if maxpw > CLUSTER_PWCUTOFF:
                    self.neighbors.append(n)

    @staticmethod
    def find_maxpw(n1, n2):
        max_val = -0.1
        for conf1 in n1.free_conformers:
            ic = conf1.i
            for conf2 in n2.free_conformers:
                jc = conf2.i
                pw = abs(pairwise[ic][jc])
                if max_val < pw:
                    max_val = pw
        return max_val


class Cluster:
    def __init__(self, nd):
        self.nodes = [nd]
        self.open_ends = False
        # Traverse to get nodes at predefined level
        queue = []
        current_node = nd
        current_level = 0
        while current_node:
            if current_level < CLUSTER_EXTENDED_LEVEL:  # only lower level needs to be extended
                for x in current_node.neighbors:
                    if x not in self.nodes:
                        queue.append((x, current_level+1))
                        self.nodes.append(x)
            else:  # at the highest level
                #print("%s:%s = %d" % (nd.resid, current_node.resid, current_level))
                for x in current_node.neighbors:
                    if x not in self.nodes:
                        self.open_ends = True     # some extended nodes are outside the level
                        break
            if self.open_ends:
               break
            if queue:
                current_node, current_level = queue.pop(0)
            else:
                self.open_ends = False    # all nodes exaused
                break
        self.accessible_states = []
        self.E_ambient = 0.0
        self.E_cluster = 0.0
        self.E_global = 0.0
        self.net_charge = 0.0
        self.dipole_direction = (0.0, 0.0, 0.0)
        self.dipole_magnitude = 0.0
        # Now get outside residues
        self.outside_residues = []
        cluster_resid = [x.resid for x in self.nodes]
        for res in residues:
            if res.resid not in cluster_resid:
                self.outside_residues.append(res)
        return

    def mc_run(self):

        return

    def update_entropy(self):
        return
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
            print("Conformer %s occurred %d times" % (name, occurrence))
            sys.exit()
    return conformers


def load_pairwise():
    folder = env.energy_table
    n_conf = len(head3lst)

    confnames = [x.confname for x in head3lst]
    scale_ele = env.var["SCALING,ELE"]
    scale_vdw = env.var["SCALING,VDW"]
    pw = np.zeros(shape=(n_conf, n_conf))
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
                pw[ic][jc] = ele * scale_ele + vdw * scale_vdw

    # average the opposite sides
    for ic in range(n_conf - 1):
        for jc in range(ic + 1, n_conf):
            # if abs(pw[ic][jc] - pw[jc][ic]) > 0.000001:
            # print("%s %.3f <-> %s %.3f" % (confnames[ic], pw[ic][jc], confnames[jc], pw[jc][ic]))
            averaged_pw = (pw[ic][jc] + pw[jc][ic]) * 0.5
            pw[ic][jc] = pw[jc][ic] = averaged_pw
    return pw


def group_residues():
    residue_ids = []
    confnames = [x.confname for x in head3lst]
    for confname in confnames:
        resid = confname[:3] + confname[5:11]
        if resid not in residue_ids:
            residue_ids.append(resid)

    residues = [Residue(x) for x in residue_ids]
    for ic in range(len(confnames)):
        confname = confnames[ic]
        resid = confname[:3] + confname[5:11]
        index = residue_ids.index(resid)
        residues[index].conformers.append(conformers[ic])

    print("   Verifying conformers ...")
    for res in residues:
        res.verify_flags()
    return residues


def report_residues():
    lines = []
    lines.append("iConf CONFORMER     FL  occ          Residue           Node           Free conformer\n")
    for res in residues:
        if len(res.free_conformers) > 0:
            iline = 0
            for conf in res.free_conformers:
                line = "%05d %s f %4.2f |" % (head3lst[conf.i].iConf, head3lst[conf.i].confname, conf.occ)
                if iline == 0:
                    line += " ----- %s ----- %s ----- | %s" % (res.resid, res.resid, head3lst[conf.i].confname)
                else:
                    line += " %s | %s" % (" " * 37, head3lst[conf.i].confname)
                lines.append("%s\n" % line)
                iline += 1
            for conf in res.fixed_conformers:
                line = "%05d %s t %4.2f |" % (head3lst[conf.i].iConf, head3lst[conf.i].confname, conf.occ)
                if iline == 0:
                    line += " ----- %s ----- %s ----- | %s" % (res.resid, res.resid, head3lst[conf.i].confname)
                lines.append("%s\n" % line)
                iline += 1
        else:
            iline = 0
            for conf in res.fixed_conformers:
                line = "%05d %s t %4.2f |" % (head3lst[conf.i].iConf, head3lst[conf.i].confname, conf.occ)
                if iline == 0:
                    line += " ----- %s" % (res.resid)
                lines.append("%s\n" % line)
                iline += 1

        lines.append("%s\n" % ("." * 84))

    open(residue_report, "w").writelines(lines)
    return


def report_neighbors():
    lines = []
    for nd in nodes:
        lines.append("%3d %s: %s\n" % (len(nd.neighbors), nd.resid , ",".join([x.resid for x in nd.neighbors])))

    lines.append("\n")
    open(neighbor_report, "w").writelines(lines)
    return

def define_nodes():
    nodes = []
    for res in residues:
        if len(res.free_conformers) > 1:
            nodes.append(Node(res))
    return nodes


def define_clusters():
    cls = []
    for nd in nodes:
        cls.append(Cluster(nd))

    # Filter clusters: merge closed clusters
    filtered_clusters = []
    while cls:
        #print("%s" % ",".join([x.nodes[0].resid for x in filtered_clusters]))
        cluster = cls.pop(0)
        filtered_clusters.append(cluster)
        if not cluster.open_ends:  # closed cluster
            # Remove any other clusters that are member of this cluster
            member_nodes = [x.resid for x in cluster.nodes[1:]]
            cls = [x for x in cls if x.nodes[0].resid not in member_nodes]
            filtered_clusters = [x for x in filtered_clusters if x.nodes[0].resid not in member_nodes]

    return filtered_clusters


def report_clusters():
    lines = []
    for cluster in clusters:
        if cluster.open_ends:
            t = "+"
        else:
            t = ""
        lines.append("%3d %s: %s %s\n" % (len(cluster.nodes), cluster.nodes[0].resid, ",".join([x.resid for x in
                                                                                        cluster.nodes]), t))

    lines.append("\n")
    open(cluster_report, "w").writelines(lines)
    return

def update_conf_energy(ph, eh):
    for ic in range(len(conformers)):
        monte_temp = env.var["(MONTE_T)"]
        E_ph = monte_temp / ROOMT * head3lst[ic].nh * (ph - head3lst[ic].pk0) * PH2KCAL
        E_eh = monte_temp / ROOMT * head3lst[ic].ne * (eh - head3lst[ic].em0) * PH2KCAL / 58.0
        conformers[ic].E_pheh = E_ph + E_eh
    return


def MC(cluster):
    """Monte Carlo sampling function."""
    # entropy run
    if env.var["(MONTE_TSX)"].upper() == "T":
        cluster.mc_run()
        cluster.update_entropy()

    # occ run
    for imc in range(6):
        cluster.mc_run()

    return



def cluster_MC(cluster):
    # initialize in-cluster conformer self energy, mfe and entropy
    for nd in cluster.nodes:
        for conf in nd.free_conformers:
            conf.E_mfe = 0.0
            for res in cluster.outside_residues:
                for conf2 in res.free_conformers:
                    conf.E_mfe += pairwise[conf.i][conf2.i] * conf2.occ
            #print("%s: %6.2f" % (head3lst[conf.i].confname, conf.E_mfe))

    # compute ambient energy
    # self
    E_self = 0.0
    for res in cluster.outside_residues:
        for conf in res.conformers:
            E_self += conf.E_self + conf.E_pheh

    # pw
    E_pw = 0.0
    n_res = len(cluster.outside_residues)
    for ires in range(n_res -1):
        for jres in range(ires+1, n_res):
            for conf1 in cluster.outside_residues[ires].conformers:
                for conf2 in cluster.outside_residues[jres].conformers:
                    E_pw += pairwise[conf1.i][conf2.i] * conf1.occ * conf2.occ

    cluster.E_ambient = E_self + E_pw

    # Test analytical
    nstate = 1
    for nd in cluster.nodes:
        nstate *= len(nd.free_conformers)

    if nstate > env.var["(NSTATE_MAX)"]:
        # Monte Carlo sampling
        if cluster.open_ends:
            t = " "
        else:
            t = "+"
        msg = "Cluster: %s %s\n" % (", ".join([x.resid for x in cluster.nodes]), t)
        fp_mc_progress.write(msg)
        msg = "n = %d, > %d Monte Carlo sampling\n" % (nstate, env.var["(NSTATE_MAX)"])
        fp_mc_progress.write(msg)
        MC(cluster)

    else:
        # Analytical solution
        if cluster.open_ends:
            t = " "
        else:
            t = "+"
        msg = "Cluster: %s %s\n" % (", ".join([x.resid for x in cluster.nodes]), t)
        fp_mc_progress.write(msg)
        msg = "n = %d, <= %d Analytical solution\n" % (nstate, env.var["(NSTATE_MAX)"])
        fp_mc_progress.write(msg)

    # Push occ to occ_old and compute new occ.
    return

def equilibrate_clusters():
    exit_eq = False
    for icycle in range(EQ_MAXCYCLES):
        q = range(len(clusters))
        random.shuffle(q)
        for icluster in q:
            cluster_MC(clusters[icluster])

        # compute convergence and energy
        # break if converged early
        break


    return



env = Env()
head3lst = load_head3lst()
conformers = [Conformer(ic) for ic in range(len(head3lst))]
pairwise = load_pairwise()
residues = group_residues()
nodes = define_nodes()
for node in nodes:
    node.find_neighbors()
clusters = define_clusters()


if __name__ == "__main__":
    report_residues()
    report_neighbors()
    report_clusters()

    ph_start = env.var["(TITR_PH0)"]
    ph_step = env.var["(TITR_PHD)"]
    eh_start = env.var["(TITR_EH0)"]
    eh_step = env.var["(TITR_EHD)"]
    titration_type = env.var["(TITR_TYPE)"]
    titration_steps = env.var["(TITR_STEPS)"]

    fp_mc_progress = open(mc_progress, "w")
    for ititr in range(titration_steps):
        if titration_type == "ph":
            ph = ph_start + ph_step * ititr
            eh = eh_start
        else:
            ph = ph_start
            eh = eh_start + eh_step * ititr

        print ("\n   Titration at pH = %.1f and Eh = %.f" % (ph, eh))
        line = "Sampling at pH = %.1f and Eh = %.f\n" % (ph, eh)
        fp_mc_progress.write(line)

        update_conf_energy(ph, eh)
        equilibrate_clusters()
        #initialize_nodes()
        #initialize_clusters(ph, eh)
    fp_mc_progress.close()
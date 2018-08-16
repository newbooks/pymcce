# pymcce

MCCE in Python




## Data structure
### Step 4: Monte Carlo sampling
env:
    (Hard define values)
    -> self.runprm = "run.prm"
    -> self.version = "PyMCCE 0.1"
    -> self.fn_conflist1 = "head1.lst"
    -> self.fn_conflist2 = "head2.lst"
    -> self.fn_conflist3 = "head3.lst"
    -> self.energy_table = "energies"

    (run.prm variables, key:value)
    -> self.var = {}

    (tpl file parameters, key1, key2, key3: value)
    -> self.param = {}



Conformer:
    ->  self.confname = fields[1]
    ->  self.flag = fields[2]
    ->  self.occ = float(fields[3])
    ->  self.crg = float(fields[4])
    ->  self.em0 = float(fields[5])
    ->  self.pk0 = float(fields[6])
    ->  self.ne = int(fields[7])
    ->  self.nh = int(fields[8])
    ->  self.vdw0 = float(fields[9]) * env.param[("SCALING", "VDW0")]
    ->  self.vdw1 = float(fields[10]) * env.param[("SCALING", "VDW1")]
    ->  self.tors = float(fields[11]) * env.param[("SCALING", "TORS")]
    ->  self.epol = float(fields[12]) * env.param[("SCALING", "ELE")]
    ->  self.dsolv = float(fields[13]) * env.param[("SCALING", "DSOLV")]
    ->  self.extra = float(fields[14])
    ->  self.history = fields[15]
    ->  self.entropy = 0.0   # -TS, will be calculated at entropy sampling
    ->  self.fixed_occ       # It starts as occ, updated with reduction


prot:
    -> self.head3list = []   (list of conformers in conformer structure, which has head3.lst terms)
    -> self.sel_energy = []  (list of self energy of each conformers, this changes with ph and entropy)
    -> self.confnames = []  (a list of conformer names for quick reference and indexing)
    -> self.pairwise = {} (pairwise interaction (i,j):ele * scaling_ele + vdw * scaling_vdw
    -> self.residues = [[]] (irregular 2d arrays that store residues in rows and conformer indices in column)
    -> self.fixed_conformers[]  (a list of fixed conformers, may have partial occ. This may vary )
    -> self.fixed_occ[] (match fixed_conformers, the occupancy of fixed conformers, using head3list and reduction)
    -> self.free_residues[[]] (free residues in rows and flippable conformers in column)
    -> self.fixed_conformers_running[]  (a copy, this data can be altered during reduced run)
    -> self.free_residues_running[[]] (a copy, this data can be altered during reduced run)
    -> self.biglist[[]] (size of free_residues_running, holds indices to res with big interaction to each residue)
    -> self.titration_type = "pH" or "Eh"
    -> self.pH  (current pH)
    -> self.Eh  (current Eh)

MicroState: (in function monte())
    -> self.state  (selection of conformers in free residues)
    -> self.fixed_conformers   (conformers that are fixed)
    -> self.complete_state   (merged state: fixed_conformers + state)
    -> self.E_state          (micro state energy)
    -> randomized_state(prot)
    -> get_E()


Monte Carlo subroutines:
    Basic function:
    mc_run(prot = prot, state = state, N=env.var["MONTE_NITER"]*len(prot.head3list), record = False)
        This will run Monte Carlo sampling N times, and record the statistics in prot.head3list[i].count

    Application:
    mc_annealing(prot = prot, start_T = env.var["MONTE_T"] + 1000, nsteps = 100)
        run mc_run() nsteps to env.var["MONTE_T"], return the state, do not record the statistics

    mc_equilibration(prot = prot, state, N = env.var["MONTE_NEQ"]*len(prot.head3list))
        starting from last state, run mc_run() with N samplings, and return the end state, do not record the statistics

    mc_reduce(prot = prot, n = 6)
        run mc_annealing(), mc_equilibration(), and mc_run() for n independent times to remove unoppcupied conformers

    mc_entropy(prot = prot, n = 6)
        run mc_annealing(), mc_equilibration(), and mc_run() for up to n independent times to calculate entropy

    mc_sample(prot = prot)
        run mc_annealing(), mc_equilibration(), and mc_run() to collect statistics

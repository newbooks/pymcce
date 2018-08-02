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

microstate: (in function monte())
    -> state[] ( a selection of conformers, one per free residue)
    -> self.complete_state[] (a merge of selected conformers composed by adding state[] and fixed_conformers[])
    -> E_state (state energy from direct calculation, sigma[complete_state * occ if partial])
    -> E_running (running state energy from differential calculation)
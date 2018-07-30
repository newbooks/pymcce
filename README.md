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

prot:
    -> self.head3list = []   (list of conformers in conformer structure, which has head3.lst terms)
    -> self.confnames = []  (a list of conformer names for quick reference and indexing)
    -> self.pairwise = {} (pairwise interaction (i,j):(ele * scaling_ele, vdw * scaling_vdw)


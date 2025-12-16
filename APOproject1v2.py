from pyomo.environ import *
import math
import matplotlib.pyplot as plt
import pandas as pd
model = ConcreteModel('APOProject')

# For Alexia to run
# inputGCs= r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx"

# For Julia to run
inputGCs= r"C:\Users\natur\Documents\AOBaddies\GCs.xlsx"

#Reading GCs from Excel sheets
df = pd.read_excel(inputGCs, sheet_name="Groups")
df.columns = df.columns.str.strip()   # Clean spaces
df['Group'] = df['Group'].astype(str).str.strip()
df = df.set_index('Group') 

dfcp=pd.read_excel(inputGCs, sheet_name="CpGC")
dfcp.columns = dfcp.columns.str.strip()   # Clean spaces
dfcp['Group'] = dfcp['Group'].astype(str).str.strip()
dfcp = dfcp.set_index('Group')  #critical for dfcp.loc[g, 'Cp']

dfSVSH = pd.read_excel(inputGCs, sheet_name="SVSH", index_col=0)
dfSVSH.columns = dfSVSH.columns.str.strip()
dfSVSH.index = dfSVSH.index.str.strip()

#Max number of times a given group can be in the molecule
#15 used in Sahinidis
ni_max=15

#Binary expension bits
Ki_bits = math.ceil(math.log2(ni_max + 1))

#----DEFINE SETS----
model.g = Set(initialize=dfcp.index.unique().tolist()) #Groups and properties for GC Cp (Sahinidis et. al)
A_props=['Tci', 'Pci', 'Tbi',  'A0i', 'B0i', 'C0i', 'D0i'] #To calculate Cp from Sahinidis
model.A=Set(initialize=A_props)

model.i = Set(initialize=df.index.unique().tolist()) #First order groups and properties for GC (Hukkerikar et. al)
X_props = ['Tb1i','Tm1i','δD1i','δP1i','δH1i','Vm1i']
model.X = Set(initialize=X_props) # The parameters calculated from MG method - might actually just be f(X) or not could put different constraints definition for some X 

model.B=Set(initialize=['acyclic','monocyclic']) #Molecule types; 'bicyclic' not included in this analysis

#Index sets for binary expansion
model.Ki = RangeSet(1, Ki_bits)

#Assigning molecule types from Excel
TYPE_COL = "Type of molecule (like aromatic and so on) ?"
type_str_dict = df[TYPE_COL].to_dict()

model.Ga = Set(initialize=[i for i, t in type_str_dict.items() if t == 'aromatic'])
print("Ga groups:", ", ".join(sorted([str(i) for i in model.Ga])))
model.Gc = Set(initialize=[i for i, t in type_str_dict.items() if t == 'cyclic'])
model.Ggen = Set(initialize=[i for i, t in type_str_dict.items() if t == 'general'])
model.Ggunsat = Set(initialize=[i for i, t in type_str_dict.items() if t == 'unsat'])

#----DEFINE PARAMETERS----
#GC for each Sahinidis group considered and each property (g,A)
g_data = {}
for g in dfcp.index:
    for A in A_props:
        g_data[(g, A)] = float(dfcp.loc[g, A])
model.cg = Param(model.g, model.A, initialize=g_data, default=0)

#GC for each Hukkerikar group considered and each property (i,X)
i_data = {}
for i in df.index:
    for X in X_props:
        i_data[(i, X)] = float(df.loc[i, X])
model.ci = Param(model.i, model.X, initialize=i_data, default=0)

#Valency data (defined for i groups)
valency_dict = df['Valency 1'].to_dict()
model.vi = Param(model.i, initialize=valency_dict, default=0)

#Mapping g groups to i groups - defined manually
ig_dict = {}
for i in dfSVSH.index:
    for g in dfSVSH.columns:
        try:
            ig_dict[(i, g)] = int(dfSVSH.loc[i, g]) if pd.notna(dfSVSH.loc[i, g]) else 0
        except:
            ig_dict[(i, g)] = 0
model.ig = Param(model.i, model.g, initialize=ig_dict, default=0)

model.Tb0 = Param(initialize=244.7889) #Kelvins
model.Tm0 = Param(initialize=144.0977) #Kelvins
model.Vm0 = Param(initialize=0.0123) #m^3/kmol
model.R0 = Param(initialize=3.3) #MPa 1/2, based on experimental work 
model.sigmacD=Param(initialize=15.7) #MPa 1/2
model.sigmacP=Param(initialize=5.2) #MPa 1/2
model.sigmacH=Param(initialize=5.8) #MPa 1/2
model.T_avg = Param(initialize=353) #in K, taken as average of absorption/desorption columns

#Scaling parameters rn using tanh but we could also use other 
#Those params are OBSOLETE
model.k_rho = Param(initialize=0.01) #Steepness for rho scaling
model.k_Cp = Param(initialize=0.05) #Steepness for Cp scaling
model.k_RED = Param(initialize=0.1) #steepness for RED scaling

model.t_rho = Param(initialize=15.828)     # MEA target for rho scaling
model.t_Cp = Param(initialize=166.6)     # MEA target for Cp scaling
model.t_RED = Param(initialize=3.94) ## MEA target for RED scaling
model.d_rho = Param(initialize=5)     # allowed change for scaling
model.d_Cp = Param(initialize=30)     # allowed change for Cp scaling
model.d_RED = Param(initialize=3) ## allowed change for Cp scaling

#Scaling params NOT USED ATM
model.rhomin=Param(initialize=0.1) 
model.rhomax=Param(initialize=36667)#in mol/m^3, estimation using literature values
model.Cpmin=Param(initialize=0.01)
model.Cpmax=Param(initialize=200) #in J/mol K, estimation using litereature values
model.REDmin=Param(initialize=0.0000001)
model.REDmax=Param(initialize=10) #in J/gK, estimation using litereature values

# Small epsilon for numerical safety
model.eps = Param(initialize=1e-12)#double check if really useful

# Objective function weights 
model.w_RED = Param(initialize=1.0)     # weight for RED
model.w_rho = Param(initialize=1.0)     # weight for density-NEGATIVE BECAUSE WE WANT TO MAXIMIZE IT
model.w_Cp  = Param(initialize=1.0)     # weight for heat capacity

#----DEFINE VARIABLES---## Need to justify bound selection
model.Tm=Var(within=NonNegativeReals, bounds=(0.01, 313), initialize=5) #in K
model.Tb=Var(within=NonNegativeReals, bounds=(393, 3000), initialize=500) #in K
model.rho=Var(within=NonNegativeReals, bounds=(0.0001, 50000)) #in mol.m^3
model.RED=Var(within=NonNegativeReals, bounds=(0.00001, None)) #in K
model.sol1=Var(within=Reals, bounds=(-16, 16)) #in 
model.sol2=Var(within=Reals, bounds=(-16, 16)) #in 
model.sol3=Var(within=Reals, bounds=(-16, 16)) #in 
model.molarvol=Var(within=NonNegativeReals, bounds=(0.000001, 100000)) #in 

#Set of binary variables (we'll later define integar cuts to get rid of some solutions, also the variables are defined through linear expressions below)
model.yi = Var(model.i, model.Ki, within=Binary)

#Binary molecule types : molecule type selection (acyclic vs monocyclic)
model.yb=Var(model.B, within=Binary) #binary molecule types defined later
model.m=Var(within=Reals, bounds=(-1,2))

#aromatic or cyclic mode 
model.ya =Var(within=Binary)  # aromatic mode Param (initialize=0) 
model.yc =Var(within=Binary)  # cyclic (non-aromatic) mode Param (initialize=0) 

#Integer counts (calculated from binary variables)
model.ni = Var(model.i, within=NonNegativeReals, bounds=(0, ni_max))#!!! watch w the integer cuts thing might have to def ni=sumk=0 to k 2^kyik (continuous associated with a set of binary variables)
model.ng = Var(model.g, within=NonNegativeReals, bounds=(0, None))

#Intermediate variables for calculating Cp, bounds copied from Sahinidis et. al where bounding doesn't result in infeasible
#Could eliminate by defining Cp but then can't bound
model.T_b = Var(within=NonNegativeReals, initialize=350, bounds=(50, 1000)) #K
model.T_c = Var(within=NonNegativeReals, initialize=500, bounds=(100, 2000)) #K
model.P_c = Var(within=PositiveReals, initialize=1.0) #bar, bounds=(2, 200)
model.T_avgr = Var(within=PositiveReals, initialize=0.7, bounds=(0.001, 1)) #K/K
model.W = Var(within=Reals, initialize=0.3) #bounds=(-1, 1.3)
model.Cp0 = Var(within=NonNegativeReals, initialize=1) #cal/mol*K (converted so Cp is in J/mol*K)

#Properties results
model.Cp = Var(within=NonNegativeReals, bounds=(1, 1000)) #J/mol*K
model.fX = Var(model.X, within=Reals)  #Tm, Tb, etc.

# Big-M for aromatic bonding
M_groups_val = ni_max * max(1, len(model.i))
model.M_groups = Param(initialize=M_groups_val)

# Indicator Params (0/1) for membership and valency > 2
is_arom_dict = {i: 1 if i in list(model.Ga) else 0 for i in model.i}
model.is_arom = Param(model.i, initialize=is_arom_dict, default=0)

is_vgt2_arom_dict = {i: 1 if (i in list(model.Ga) and int(valency_dict.get(i, 0)) > 2) else 0 for i in model.i}
model.is_vgt2_arom = Param(model.i, initialize=is_vgt2_arom_dict, default=0)

# Expressions for counts
model.count_arom_vgt2 = Expression(expr=sum(model.ni[i] * model.is_vgt2_arom[i] for i in model.i))
model.count_non_aromatic = Expression(expr=sum(model.ni[i] * (1 - model.is_arom[i]) for i in model.i))

# Binary indicating whether attachment out of the aromatic ring is permitted
model.z_attach_ok = Var(within=Binary)

#---DEFINE CONSTRAINTS---#
# (1) Property group contribution relationship
def XGC_rule(model, X):
    return model.fX[X]==sum(model.ni[i] * model.ci[i, X] for i in model.i)
model.XGC=Constraint(model.X,rule=XGC_rule)

# (2) Number of g-groups derived from i-groups
def numgini_rule(model,g):#either we define specific exceptions to take into account the various f(X) or idk
    return model.ng[g]==sum(model.ni[i] * model.ig[i,g] for i in model.i)
model.numgini = Constraint(model.g, rule=numgini_rule)

# (3) Cp calculations from Sahinidis et al 
# Critical temperature/pressure, boiling temperature from set of g-groups in final molecule
# Boiling temperature
def Tb_calc_rule(model):
    return model.T_b == 198.2 + sum(model.ng[g] * model.cg[g, 'Tbi'] for g in model.g)
model.Tb_calc = Constraint(rule=Tb_calc_rule)

# Critical temperature (reformulated to avoid division)
def Tc_calc_rule(model):
    return model.T_c * (0.584 + 0.965 * sum(model.ng[g] * model.cg[g, 'Tci'] for g in model.g)) == model.T_b
model.Tc_calc = Constraint(rule=Tc_calc_rule)

# Critical pressure
def Pc_calc_rule(model):
    return model.P_c == (0.113 + 0.0032 * sum(model.ng[g] for g in model.g) - sum(model.ng[g] * model.cg[g, 'Pci'] for g in model.g))
model.Pc_calc = Constraint(rule=Pc_calc_rule)

# Reduced average temperature (reformulated)
def Tavgr_calc_rule(model):
    return model.T_avgr * model.T_c == model.T_avg
model.Tavgr_calc = Constraint(rule=Tavgr_calc_rule)

# Acentric factor
def W_calc_rule(model):
    alpha = (-5.97214 - log(model.P_c / 1.013 + (6.09648 * model.T_c) / model.T_b) + 
             1.28862 * log(model.T_b / model.T_c) - 0.167347 * (model.T_b / model.T_c)**6)
    beta = (15.2518 - (15.6875 * model.T_c) / model.T_b - 13.4721 * log(model.T_b / model.T_c) + 
            0.43577 * (model.T_b / model.T_c)**6)
    return model.W * beta == alpha
model.W_calc = Constraint(rule=W_calc_rule)

# Ideal liquid heat capacity
def Cp0_calc_rule(model):
    return model.Cp0 == (sum(model.ng[g] * model.cg[g, 'A0i'] for g in model.g) - 37.0/93.0 +
                         (sum(model.ng[g] * model.cg[g, 'B0i'] for g in model.g) + 0.21) * model.T_avg +
                         (sum(model.ng[g] * model.cg[g, 'C0i'] for g in model.g) - 3.91e-4) * model.T_avg**2 +
                         (sum(model.ng[g] * model.cg[g, 'D0i'] for g in model.g) + 2.06e-7) * model.T_avg**3)
model.Cp0_calc = Constraint(rule=Cp0_calc_rule)

# Final Cp liquid calculation (Rowlinson)
def Cp_rule(model):#either we define specific exceptions to take into account the various f(X) or idk
    return model.Cp == 4.184 * ((1/4.1868)*(model.Cp0 + 8.314*(1.45 + (0.45/(1-model.T_avgr)) + 0.25*model.W * (17.11 + 25.2*(((1 - model.T_avgr)**(1/3))/model.T_avgr) + (1.742 / (1-model.T_avgr))))))
model.Cp_constraint = Constraint(rule=Cp_rule)

# (4) Def the various functions in Hurekkikar group contributions, not sure it properly relates to our properties definitions
def link_properties_rule(model, X):
    if X == 'Tm1i':
        print(model.fX[X])
        return model.fX[X] == exp(model.Tm/model.Tm0)
    elif X == 'Tb1i':
        return model.fX[X] == exp(model.Tb/model.Tb0)
    elif X == 'Vm1i':
        return model.molarvol == model.fX[X] + model.Vm0
    elif X == 'δD1i':
        return model.sol1 == model.fX[X]
    elif X == 'δP1i':
        return model.sol2 == model.fX[X]
    elif X == 'δH1i':
        return model.sol3 == model.fX[X]
    else:
        return Constraint.Skip
model.link_properties = Constraint(model.X, rule=link_properties_rule)

# (5) Definition of density from GCS calculated parameters
def rho_def_rule(model):
    return model.rho * model.molarvol == 1 
model.rho_def = Constraint(rule=rho_def_rule)

#(6) Def of RED from GCS calc parameters
def RED_def_rule(model):
    return model.RED * model.R0 == sqrt(
        4 * (model.sol1 - model.sigmacD)**2 +
        (model.sol2 - model.sigmacP)**2 +
        (model.sol3 - model.sigmacH)**2
    )
model.RED_def = Constraint(rule=RED_def_rule)

#(7) Allowing only one type of molecule
def one_type(model):
    return sum (model.yb[b] for b in model.B) ==1
model.one_type = Constraint(rule=one_type)

#(8) Relating m to binary molecule types:
def m_rule(model):
    return model.m == 1-model.yb['monocyclic']
model.m_rule = Constraint(rule=m_rule)

#(9) Replacing integer variables with linear expressions: i 
def ni_rule(model,i):
    return model.ni[i] == sum(2**(k-1)*model.yi[i,k] for k in model.Ki)
model.ni_rule = Constraint(model.i,rule=ni_rule)

#(9.5_integer cut) to exclude a specific composition 
# Target counts to exclude:
excl_combo = {
    'CH3': 1,
    'NH2 except as above': 1,
    '-O-': 2,
}

# Utility: convert an integer count into binary bits aligned with your weights 2^(k-1)
def _bits_for_count(count, Ki_bits):
    bits = {}
    remaining = int(count)
    for k in range(1, Ki_bits + 1):
        w = 2**(k-1)
        if remaining >= w:
            bits[k] = 1
            remaining -= w
        else:
            bits[k] = 0
    return bits

# Build required bit pattern per group (fail fast if name not found)
Ki_bits_val = int(model.Ki.last())  # safer than max(model.Ki)
excl_bits = {}
for i_name, c in excl_combo.items():
    if i_name not in list(model.i):
        raise ValueError(f"Nogood cut: group '{i_name}' not found in model.i. Available: {list(model.i)}")
    excl_bits[i_name] = _bits_for_count(c, Ki_bits_val)

# Define the nogood constraint: sum of mismatches across all specified bits >= 1
def intcut_rule(model):
    return sum(
        (1 - model.yi[i, k]) if excl_bits[i][k] == 1 else model.yi[i, k]
        for i in excl_bits.keys()
        for k in model.Ki
    ) >= 1

model.intcut = Constraint(rule=intcut_rule)

'''
#(10) Replacing integer variables with linear expressions : g ###Maybe obsolete ?
def ng_rule(model,g):
    return model.ng[g] == sum(2**(k-1)*model.yg[g,k] for k in model.Kg)
model.ng_rule = Constraint(model.g, rule=ng_rule)
'''
#(11) Octet Rule, i compounds
def vi_rule(model):
    return 2*model.m==sum((2-model.vi[i])*model.ni[i] for i in model.i)
model.vi_rule = Constraint(rule=vi_rule)

#(12) One cyclic mode
# Exactly one mode if monocyclic, none if acyclic
def monocyclic_mode_selection_rule(model):
    return model.ya + model.yc == model.yb['monocyclic']
model.monocyclic_mode_selection = Constraint(rule=monocyclic_mode_selection_rule)

#(13) Aromaticity constratints
# Exactly 6 aromatic groups if aromatic mode; 0 aromatic groups otherwise
def aromatic_ring_rule(model):
    return sum(model.ni[i] for i in model.Ga) == 6 * model.ya
model.aromatic_ring = Constraint(rule=aromatic_ring_rule)

# Forbid cyclic groups unless cyclic mode is chosen (Big-M)
def cyclic_only_if_cyclic_mode_rule(model):
    return sum(model.ni[i] for i in model.Gc) <= model.M_groups * model.yc
model.cyclic_only_if_cyclic_mode = Constraint(rule=cyclic_only_if_cyclic_mode_rule)

#(14) Aromatic group with a valency of >2 only can bind to non-aromatic groups
#count expression
def attach_ok_lower_rule(model):
    return model.count_arom_vgt2 >= model.z_attach_ok
model.attach_ok_lower = Constraint(rule=attach_ok_lower_rule)

def attach_ok_upper_rule(model):
    return model.count_arom_vgt2 <= model.M_groups * model.z_attach_ok
model.attach_ok_upper = Constraint(rule=attach_ok_upper_rule)

def non_aromatic_allowed_in_aromatic_mode_rule(model):
    return model.count_non_aromatic <= model.M_groups * (1 - model.ya) + model.M_groups * model.z_attach_ok
model.non_aromatic_allowed_in_aromatic_mode = Constraint(rule=non_aromatic_allowed_in_aromatic_mode_rule)

#(15) Bonding rule 
def bi_rule(model, i):
    return model.ni[i]*(model.vi[i]-1)+2-sum(model.ni[u] for u in model.i)<=0 
model.bi_rule = Constraint(model.i, rule=bi_rule)

#(16) Minimum of two groups bonded to each other
def zer(model):
    return sum(model.ni[k] for k in model.i)>=2
model.zer_rule = Constraint(rule=zer)    

#(17) Maximum amount of molecules
def maxgroups(model):
    return sum(model.ni[k] for k in model.i)<=25
model.max_rule = Constraint(rule=maxgroups)    

#(18) Minimum of 5 cyclic groups if cyclic non-aromatic mode
def min_cyclic_groups_rule(model):
    return sum(model.ni[i] for i in model.Gc) >= 5 * model.yc
model.min_cyclic_groups = Constraint(rule=min_cyclic_groups_rule)

#---DEFINE OBJECTIVE FUNCTION---#
'''
#Normalized Expression using traditional scaling
model.rho_norm = Expression(
    expr=(model.rho - model.rhomin) / (model.rhomax - model.rhomin))
model.Cp_norm = Expression(
    expr=(model.Cp - model.Cpmin) / (model.Cpmax - model.Cpmin))
model.RED_norm = Expression(
    expr=(log(model.RED) - log(model.REDmin)) / (log(model.REDmax) - log(model.REDmin)))
'''
#standardized data using target
model.rho_norm = Expression(
    expr=(model.rho - model.t_rho) / (model.d_rho))
model.Cp_norm = Expression(
    expr=(model.Cp - model.t_Cp) / (model.d_Cp))
model.RED_norm = Expression(
    expr=(model.RED - model.t_RED) / (model.d_RED))

#Actual objective function
def objective_rule(model):
    return (
        model.w_RED * model.RED_norm +
        model.w_rho *(1- model.rho_norm) +
        model.w_Cp  * model.Cp_norm
    )
model.obj = Objective(rule=objective_rule, sense=minimize)

#----SOLVE----#
solver = SolverFactory('gams')
solver.options['solver'] = 'DICOPT'
results = solver.solve(model, tee=True)

#----DISPLAY RESULTS----#
print("\n" + "="*70)
print("OPTIMAL SOLUTION FOR CARBON CAPTURE SOLVENT DESIGN")
print("="*70)

print("\nGroup Counts (ni) - Hukkerikar Groups:")
print("-" * 70)
for i in model.i:
    if model.ni[i].value and model.ni[i].value > 0.01:
        print(f"  {i:25s}: {model.ni[i].value:6.3f}")

print("\nSecondary Group Counts (ng) - Sahinidis Groups:")
print("-" * 70)
for g in model.g:
    if model.ng[g].value and model.ng[g].value > 0.01:
        print(f"  {g:25s}: {model.ng[g].value:6.3f}")

print("\nMolecule Type:")
print("-" * 70)
for b in model.B:
    if model.yb[b].value and model.yb[b].value > 0.5:
        print(f"  Type: {b}")
        print(f"  m parameter: {model.m.value:.3f}")

if model.ya.value and model.ya.value > 0.5:
    print("  Cyclic Mode: Aromatic")
if model.yc.value and model.yc.value > 0.5:
    print("  Cyclic Mode: Non-aromatic cyclic")

print("\nTarget Properties:")
print("-" * 70)
print(f"  Heat Capacity (Cp):        {model.Cp.value:.6f} J/(mol·K)")
print(f"  Density (rho):             {model.rho.value:.6f} mol/m³")
print(f"  RED (Solubility):          {model.RED.value:.6f}")
print(f"  Melting Point (Tm):        {model.Tm.value:.3f} K ({model.Tm.value - 273.15:.2f} °C)")
print(f"  Boiling Point (Tb):        {model.Tb.value:.3f} K ({model.Tb.value - 273.15:.2f} °C)")
print(f"  Molar Volume:              {model.molarvol.value:.4f} m³/kmol")

print("\nIntermediate Cp Calculation Values:")
print("-" * 70)
if model.T_b.value:
    print(f"  Boiling Temperature (T_b): {model.T_b.value:.2f} K")
if model.T_c.value:
    print(f"  Critical Temperature:      {model.T_c.value:.2f} K")
if model.P_c.value:
    print(f"  Critical Pressure:         {model.P_c.value:.4f} bar")
if model.T_avgr.value:
    print(f"  Reduced Avg Temperature:   {model.T_avgr.value:.4f}")
if model.W.value:
    print(f"  Acentric Factor (W):       {model.W.value:.4f}")
if model.Cp0.value:
    print(f"  Ideal Cp0:                 {model.Cp0.value:.6f} cal/(mol·K)")

print("\nSolubility Parameters:")
print("-" * 70)
print(f"  δD (Dispersion):           {model.sol1.value:.4f} MPa^0.5")
print(f"  δP (Polar):                {model.sol2.value:.4f} MPa^0.5")
print(f"  δH (Hydrogen bonding):     {model.sol3.value:.4f} MPa^0.5")

print("\nObjective Function Components:")
print("-" * 70)
print(f"  RED normalized:            {model.RED_norm():.6f}")
print(f"  Density normalized:        {model.rho_norm():.6f}")
print(f"  Cp normalized:             {model.Cp_norm():.6f}")
print(f"  Total Objective Value:     {model.obj():.6f}")

print("\n" + "="*70)
print("Solver Status:", results.solver.status)
print("Termination Condition:", results.solver.termination_condition)
print("="*70)
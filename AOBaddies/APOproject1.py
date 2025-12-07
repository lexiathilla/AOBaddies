from pyomo.environ import *
import math
import matplotlib.pyplot as plt
import pandas as pd
model = ConcreteModel('APOProject')
inputGCs="C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx"

#Excel Reads
df = pd.read_excel(inputGCs, sheet_name="Group")
df.columns = df.columns.str.strip()   # Clean spaces

dfcp=pd.read_excel(inputGCs, sheet_name="CpGC")
dfcp.columns = dfcp.columns.str.strip()   # Clean spaces

dfSVSH = pd.read_excel("GCs.xlsx", sheet_name="SVSH", index_col=0)
dfSVSH.columns = dfSVSH.columns.str.strip()
dfSVSH.index = dfSVSH.index.str.strip()

#max number of groups per elements
ni_max=8
ng_max=10

#Define SETS
model.g = Set(initialize=dfcp['Group'].unique().tolist())#groups for GC Cp
model.i = Set(initialize=df['Group'].unique().tolist())#First order groups for GC - need to choose from large variety - litterature backed ?
model.A=Set(initialize=['Cp'])
X_props = ['Tb1i','Tm1i','δD1i','δP1i','δH1i','Vm1i']
model.X = Set(initialize=X_props)#The parameters calculated from MG method - might actually just be f(X) or not could put different constraints definition for some X 
model.B=Set(initialize=['acyclic','monocyclic'])#No'bicyclic'

#Ki={i: math.ceil(math.log2(ni)) for i in model.i}  could create ngs and nis for each g and i btw
#Kg={g: math.ceil(math.log2(ng)) for g in model.g}

# Bit-lengths for binary expansion (correct Set initialization)
Ki = math.ceil(math.log2(ni_max))
Kg = math.ceil(math.log2(ng_max))
# Create index sets
model.Ki = Set(1, Ki + 1)
model.Kg = Set(1, Kg + 1)

#Define Parameters GC, ideally for each grp considered and each property (i,X)
i_data = {}
for _, row in df.iterrows():
    i = row['Group']
    for X in X_props:
        i_data[(i, X)] = float(row[X])
model.ci = Param(model.i, model.X, initialize=i_data, default=0)

#for Cp, group contributions
g_data = { (g, 'Cp'): float(dfcp.loc[g, 'Cp'])
           for g in model.g if 'Cp' in dfcp.columns }
model.cg = Param(model.g, model.A, initialize=g_data, default=0)


#valency data (defined for i for now)
valency_dict = df.set_index('Group')['Valency 1'].to_dict()
model.vi = Param(model.i, initialize=valency_dict, default=0)
#v_g={}
#model.vg=Param(model.i, initialize=v_g, default=0)

#group types mapping ; probably wrong way to do it : Just directly use column
TYPE_COL = "Type of molecule (like aromatic and so on) ?"

# Build a dict: Group -> type string (read exactly as in Excel)
type_str_dict = df.set_index('Group')[TYPE_COL].to_dict()

# Create Pyomo sets from the dict (no normalization)
model.Ga = Set(initialize=[i for i, t in type_str_dict.items() if t == 'aromatic'])  # aromatic groups
model.Gc = Set(initialize=[i for i, t in type_str_dict.items() if t == 'cyclic'])    # cyclic groups
model.Ggen = Set(initialize=[i for i, t in type_str_dict.items() if t == 'general']) # general groups
model.Ggunsat = Set(initialize=[i for i, t in type_str_dict.items() if t == 'unsat']) # general groups
#Number of groups in other groups (link between Cp and Hurekkikar et al) #TO CORRECT, when relation not mentionned we get zero
ig_dict = {
    (i, g): int(dfSVSH.loc[i, g])
    for i in dfSVSH.index
    for g in dfSVSH.columns
}
model.ig = Param(model.i, model.g, initialize=ig_dict, default=0)

#----DEFINE VARIABLES--- ##WORK ON PROPERLY DEFINING BOUNDS
model.Tm=Var(within=NonNegativeReals, bounds=(0.01, 313)) #in K
model.Tb=Var(within=NonNegativeReals, bounds=(0.01, 393)) #in K
model.rho=Var(within=NonNegativeReals, bounds=(0.0001, None)) #in DEFINE
model.RED=Var(within=NonNegativeReals, bounds=(0.001, None)) #in K
model.sol1=Var(within=NonNegativeReals, bounds=(0.0001, None)) #in K
model.sol2=Var(within=NonNegativeReals, bounds=(0.0001, None)) #in K
model.sol3=Var(within=NonNegativeReals, bounds=(0.0001, None)) #in K
model.molarvol=Var(within=NonNegativeReals, bounds=(0.0001, None)) #in K

#Set of binary variables (we'll later define integar cuts to get rid of some solutions, also the variables are defined through linear expressions below)
model.yi = Var(model.i, model.Ki, within=Binary)
model.yg = Var(model.g, model.Kg, within=Binary)

#Binary molecule types : molecule type selection (acyclic vs monocyclic)
model.yb=Var(model.B, within=Binary) #binary molecule types defined later
model.m=Var(within=Reals)

#aromatic or cyclic mode 
model.ya = Var(within=Binary)  # aromatic mode
model.yc = Var(within=Binary)  # cyclic (non-aromatic) mode

#Integer counts (calculated with binary variables)
model.ni=Var(model.i, within=NonNegativeReals)#!!! watch w the integer cuts thing might have to def ni=sumk=0 to k 2^kyik (continuous associated with a set of binary variables)
model.ng = Var(model.g, within=NonNegativeReals)

model.Cp = Var(within=NonNegativeReals)
model.fX = Var(model.X, within=NonNegativeReals)  # property results (Tm, Tb, etc.)

#---PARAMETERS---
model.Tb0 = Param(initialize=244.7889)#Kelvins
model.Tm0 = Param(initialize=144.0977)#Kelvins
model.Vm0 = Param(initialize=20.7339)#m^3/kmol
model.R0 = Param(initialize=4.7)#MPa 1/2, based on experimental work 
model.sigmacD=Param(initialize=15.6)#MPa 1/2
model.sigmacP=Param(initialize=5.2)#MPa 1/2
model.sigmacH=Param(initialize=5.8)#MPa 1/2

## --- scaling parameters --- rn using tanh but we could also use other
model.k_rho = Param(initialize=0.01)     # steepness for rho scaling
model.k_Cp = Param(initialize=0.05)     # steepness for Cp scaling
model.k_RED = Param(initialize=0.1) ## steepness for RED scaling

# --- Objective function weights ---
model.w_RED = Param(initialize=1.0)     # weight for RED
model.w_rho = Param(initialize=-1.0)     # weight for density-NEGATIVE BECAUSE WE WANT TO MAXIMIZE IT
model.w_Cp  = Param(initialize=1.0)     # weight for heat capacity

# --- Constraints ---
# (1) property group contribution relationship
def XGC_rule(model, X):
    return model.fX[X]==(sum(model.ni[i] * model.ci[i, X]) for i in model.i)
model.XGC=Constraint(model.X,rule=XGC_rule)

# (2) Number of g-groups derived from i-groups
def numgini_rule(model,g):#either we define specific exceptions to take into account the various f(X) or idk
    return model.ng[g]==sum(model.ni[i] * model.ig[i,g] for i in model.i)
model.numgini = Constraint(model.g, rule=numgini_rule)

# (3) Cp calculations from g-groups
def CpGC_rule(model):#either we define specific exceptions to take into account the various f(X) or idk
    return model.Cp == sum(model.ng[g] * model.cg[g, 'Cp'] for g in model.g)
model.CpGC = Constraint(rule=CpGC_rule)

# (4) def the various functions in Hurekkikar group contributions, not sure it properly relates to our properties definitions
def link_properties_rule(model, X):
    if X == 'Tm1i':
        return model.Tm == exp(model.fX[X] / model.Tm0)
    elif X == 'Tb1i':
        return model.Tb == exp(model.fX[X] / model.Tb0)
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
    return model.rho == 1 / model.molarvol
model.rho_def = Constraint(rule=rho_def_rule)

#(6) Def of RED from GCS calc parameters
def RED_def_rule(model):
    return model.RED == sqrt(
        4 * (model.sol1 - model.sigmacD)**2 +
        (model.sol2 - model.sigmacP)**2 +
        (model.sol3 - model.sigmacH)**2
    ) / model.R0
model.RED_def = Constraint(rule=RED_def_rule)

#(7) Allowing only one type of molecule
def one_type(model):
    return sum (model.yb[b] for b in model.B) ==1
model.one_type = Constraint(rule=one_type)

#(8) Relating m to binary molecule types:
def m_rule(model):
    return model.m == model.yb['monocyclic']
model.m_rule = Constraint(rule=m_rule)

#(9) Replacing integer variables with linear expressions : i 
def ni_rule(model,i):
    return model.ni[i] == sum(2**(k-1)*model.yi[i,k] for k in model.Ki)
model.ni_rule = Constraint(rule=ni_rule)

#(10) Replacing integer variables with linear expressions : g ###Maybe obsolete ?
def ng_rule(model,g):
    return model.ng[g] == sum(2**(k-1)*model.yg[g,k] for k in model.Kg)
model.ng_rule = Constraint(rule=ng_rule)

#(11) Octet Rule, i compounds
def vi_rule(model):
    return 2*model.m==sum((2-model.vi[i])*model.ni[i] for i in model.i)
model.vi_rule = Constraint(rule=vi_rule)

#(12) one cyclic mode
# Exactly one mode if monocyclic, none if acyclic
def monocyclic_mode_selection_rule(model):
    return model.ya + model.yc == model.yb['monocyclic']
model.monocyclic_mode_selection = Constraint(rule=monocyclic_mode_selection_rule)

#(13)aromaticity
# Exactly 6 aromatic groups if aromatic mode; 0 aromatic groups otherwise
def aromatic_ring_rule(model):
    return sum(model.ni[i] for i in model.Ga) == 6 * model.ya
model.aromatic_ring = Constraint(rule=aromatic_ring_rule)

# Forbid cyclic groups unless cyclic mode is chosen (Big-M)
M_groups = ni_max * max(1, len(model.i))  # safe upper bound
model.M_groups = Param(initialize=M_groups)

def cyclic_only_if_cyclic_mode_rule(model):
    return sum(model.ni[i] for i in model.Gc) <= model.M_groups * model.yc
model.cyclic_only_if_cyclic_mode = Constraint(rule=cyclic_only_if_cyclic_mode_rule)

#Normalized expressions using tanh(kx)
model.rho_norm = Expression(expr=tanh( model.rho * model.k_rho))##maybe not we want the larger it is the smaller it is so could put a sorth of inverse there instead of w the weights
model.Cp_norm  = Expression(expr=tanh( model.Cp * model.k_Cp))
model.RED_norm  = Expression(expr=tanh(model.RED * model.k_RED))

##Actual objective function
def objective_rule(model):
    return (
        model.w_RED * model.RED_norm +
        model.w_rho * model.rho_norm +
        model.w_Cp  * model.Cp_norm
    )
model.obj = Objective(rule=objective_rule, sense=minimize)
#Solve
solver=SolverFactory('gams:ANTIGONE')#NOT sure if properly defined
results=solver.solve(model, tee=True)

#Display results
# --- Display results ---
print("\n--- Optimal Solution ---")

# --- Optional: plot ---
print("\n--- Optimal Solution ---")
for i in model.i:
    print(f"{i}: ni = {model.ni[i].value}")
print(f"Cp = {model.Cp.value}")
print(f"rho = {model.rho.value}")
print(f"RED = {model.RED.value}")


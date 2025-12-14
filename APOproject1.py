from pyomo.environ import *
import math
import matplotlib.pyplot as plt
import pandas as pd
model = ConcreteModel('APOProject')
inputGCs= r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx"

#Excel Reads
df = pd.read_excel(inputGCs, sheet_name="Group")
df.columns = df.columns.str.strip()   # Clean spaces

dfcp=pd.read_excel(inputGCs, sheet_name="CpGC")
dfcp.columns = dfcp.columns.str.strip()   # Clean spaces
dfcp['Group'] = dfcp['Group'].astype(str).strip()  # ensure clean labels
dfcp = dfcp.set_index('Group') 

dfSVSH = pd.read_excel(inputGCs, sheet_name="SVSH", index_col=0)
dfSVSH.columns = dfSVSH.columns.str.strip()
dfSVSH.index = dfSVSH.index.str.strip()

#max number of groups per elements
ni_max=8
ng_max=10

#Define SETS
model.g = Set(initialize=dfcp.index.unique().tolist())#groups for GC Cp
model.i = Set(initialize=df['Group'].unique().tolist())#First order groups for GC - need to choose from large variety - litterature backed ?
model.A=Set(initialize=['Cp'])
X_props = ['Tb1i','Tm1i','δD1i','δP1i','δH1i','Vm1i']
model.X = Set(initialize=X_props)#The parameters calculated from MG method - might actually just be f(X) or not could put different constraints definition for some X 
model.B=Set(initialize=['acyclic','monocyclic'])#No'bicyclic'

#Ki={i: math.ceil(math.log2(ni)) for i in model.i}  could create ngs and nis for each g and i btw
#Kg={g: math.ceil(math.log2(ng)) for g in model.g}

# Bit-lengths for binary expansion (correct Set initialization)
Ki_bits = math.ceil(math.log2(ni_max + 1))
Kg_bits = math.ceil(math.log2(ng_max + 1))
model.Ki = RangeSet(1, Ki_bits)
model.Kg = RangeSet(1, Kg_bits)

#Define Parameters GC, ideally for each grp considered and each property (i,X)
i_data = {}
for _, row in df.iterrows():
    i = row['Group']
    for X in X_props:
        i_data[(i, X)] = float(row[X])
model.ci = Param(model.i, model.X, initialize=i_data, default=0)

#for Cp, group contributions
#g_data = { (g, 'Cp'): float(dfcp.loc[g, 'Cp'])
#           for g in model.g if 'Cp' in dfcp.columns }
#model.cg = Param(model.g, model.A, initialize=g_data, default=0)

g_map = dfcp['Cp'].astype(float).to_dict()
g_data = {(g, 'Cp'): g_map.get(g, 0.0) for g in model.g}

model.cg = Param(model.g, model.A, initialize=g_data, default=0.0)

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
model.RED=Var(within=NonNegativeReals, bounds=(0.0001, None)) #in K
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
M_groups_vals = ni_max * max(1, len(model.i))
model.M_groups = Param(initialize=M_groups_vals)
model.ng = Var(model.g, within=NonNegativeReals)

model.Cp = Var(within=NonNegativeReals)
model.fX = Var(model.X, within=NonNegativeReals)  # property results (Tm, Tb, etc.)

#Handling aromatic bonding
# Big-M for counts (upper bound on total non-aromatic groups)
if not hasattr(model, 'M_groups'):


# Indicator Params (0/1) for membership and valency > 2
is_arom_dict = {i: 1 if i in list(model.Ga) else 0 for i in model.i}
is_vgt2_arom_dict = {i: 1 if (i in list(model.Ga) and int(valency_dict.get(i, 0)) > 2) else 0 for i in model.i}

model.is_arom = Param(model.i, initialize=is_arom_dict, default=0)
model.is_vgt2_arom = Param(model.i, initialize=is_vgt2_arom_dict, default=0)

# Expressions for counts
model.count_arom_vgt2 = Expression(expr=sum(model.ni[i] * model.is_vgt2_arom[i] for i in model.i))
model.count_non_aromatic = Expression(expr=sum(model.ni[i] * (1 - model.is_arom[i]) for i in model.i))

# Binary indicating whether attachment out of the aromatic ring is permitted
model.z_attach_ok = Var(within=Binary)

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
model.rhomin=Param(initialize=9000)
model.rhomax=Param(initialize=36667)#in mol/m^3 estimation using litterature values
model.Cpmin=Param(initialize=0.01)
model.Cpmax=Param(initialize=0.167)#in J/mol K, estimation using littereature values as well
model.REDmin=Param(initialize=0.0001)
model.REDmax=Param(initialize=1)#in J/gK, estimation using littereature values as well
# small epsilon for numerical safety
model.eps = Param(initialize=1e-12)

# --- Objective function weights ---
model.w_RED = Param(initialize=1.0)     # weight for RED
model.w_rho = Param(initialize=1.0)     # weight for density-NEGATIVE BECAUSE WE WANT TO MAXIMIZE IT
model.w_Cp  = Param(initialize=1.0)     # weight for heat capacity

# --- Constraints ---
# (1) property group contribution relationship
def XGC_rule(model, X):
    return model.fX[X]==sum(model.ni[i] * model.ci[i, X] for i in model.i)
model.XGC=Constraint(model.X,rule=XGC_rule)

# (2) Number of g-groups derived from i-groups
def numgini_rule(model,g):#either we define specific exceptions to take into account the various f(X) or idk
    return model.ng[g]==sum(model.ni[i] * model.ig[i,g] for i in model.i)
model.numgini = Constraint(model.g, rule=numgini_rule)

# (3) Cp calculations from Sahinidis et al
# # TO CORRECT - need to fix to be model.var, and to index through corresponding rows of model.cg

# Critical temperature/pressure, boiling temperature from set of g-groups in final molecule
T_b = 198.2 + sum(model.ng[i]*Tbi[i])
T_c = T_b/(0.584 + 0.965*sum(model.ng[i]*model.cg['Tci']))
P_c = (0.113 + 0.0032*sum(model.ng[i]) - sum(model.ng[i]*model.cg['Pci']))

# Reduced average temperature
T_avgr = model.T_avg/T_c

# Acentric factor
alpha = -5.97214 - log((P_c / 1.013) + ((6.09648 * T_c) / T_b)) + 1.28862*log(T_b/T_c) - 0.167347*(T_b/T_c)^6
beta = 15.2518 - ((15.6875*T_c) / T_b) - 13.4721*log(T_b/T_c) + 0.43577*(T_b / T_c)^6
model.W = alpha/ / beta

# Ideal liquid heat capacity -- WHICH NEED TO BE DEFINED AS RULES??
model.Cp0 == sum(model.ng[i]*A0i[i]) - 37/93 + (sum(model.ng[i]*B0i[i])+0.21)*model.T_avg + (sum(model.ng[i]*C0i[i])-3.91e-4)*model.T_avg^2 + (sum(model.ng[i]*D0i[i])+2.06e-7)*model.T_avg^3

# final Cp liquid calculation (Rowlinson)
def Cp_rule(model):#either we define specific exceptions to take into account the various f(X) or idk
    return model.Cp == (1/4.1868)*(model.Cp0 + 8.314(1.45 + (0.45/(1-T_avgr)) + 0.25*model.W * (17.11 + 25.2 (((1 - T_avgr)^(1/3))/T_avgr) + (1.742 / (1-T_avgr)))))
model.Cp_con= Constraint(rule=Cp_rule)


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
def cyclic_only_if_cyclic_mode_rule(model):
    return sum(model.ni[i] for i in model.Gc) <= model.M_groups * model.yc
model.cyclic_only_if_cyclic_mode = Constraint(rule=cyclic_only_if_cyclic_mode_rule)

#(14) aromatic group with a valency of >2 only can bind to non-aromatic groups

# Link z_attach_ok to presence of at least one aromatic group with valency > 2
def attach_ok_lower_rule(model):
    return model.count_arom_vgt2 >= model.z_attach_ok
model.attach_ok_lower = Constraint(rule=attach_ok_lower_rule)

def attach_ok_upper_rule(model):
    return model.count_arom_vgt2 <= model.M_groups * model.z_attach_ok
model.attach_ok_upper = Constraint(rule=attach_ok_upper_rule)

# Gate non-aromatic groups by aromatic mode and attach_ok
def non_aromatic_allowed_in_aromatic_mode_rule(model):
    return model.count_non_aromatic <= model.M_groups * (1 - model.ya) + model.M_groups * model.z_attach_ok
model.non_aromatic_allowed_in_aromatic_mode = Constraint(rule=non_aromatic_allowed_in_aromatic_mode_rule)

#Normalized expressions using tanh(kx)
#model.rho_norm = Expression(expr=tanh( model.rho * model.k_rho))##maybe not we want the larger it is the smaller it is so could put a sorth of inverse there instead of w the weights
#model.Cp_norm  = Expression(expr=tanh( model.Cp * model.k_Cp))
#model.RED_norm  = Expression(expr=tanh(model.RED * model.k_RED))#why is RED white here and not the other guys
#Normalized Expression using traditional scaling
model.rho_norm = Expression(
    expr=(model.rho - model.rhomin) / (model.rhomax - model.rhomin + model.eps))
model.Cp_norm = Expression(
    expr=(model.Cp - model.Cpmin) / (model.Cpmax - model.Cpmin + model.eps))
model.RED_norm = Expression(
    expr=(log(model.RED) - log(model.REDmin)) / (log(model.REDmax) - log(model.REDmin) + model.eps))
##Actual objective function
def objective_rule(model):
    return (
        model.w_RED * model.RED_norm +
        model.w_rho *(1- model.rho_norm) +
        model.w_Cp  * model.Cp_norm
    )
model.obj = Objective(rule=objective_rule, sense=minimize)

#Solve
solver = SolverFactory('gams')
solver.options['solver'] = 'antigone'
results = solver.solve(model, tee=True)

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


from pyomo.environ import *
from pyomo.opt import TerminationCondition, SolverStatus, check_optimal_termination
import math
import matplotlib.pyplot as plt
import pandas as pd
import os

model = ConcreteModel('APOProject')##DO WE NEED THIS LINE : I don't think we do ?

#USER OPTIONS
specni=True #set excel determined target nis (True) or default (False)
ni_max_default=15 # default if specni=False or data missing 15 used in Sahinidis et. al
most=25 #max total groups
least=2 #min total groups

# For Alexia to run
inputGCs= r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx"

# For Julia to run
#inputGCs= r"C:\Users\natur\Documents\AOBaddies\GCs.xlsx"

#READING DATA
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

#IF NEEDED, making sure ni is in the correct format
'''
# Ensure ni limits fields exist or sensible defaults
ni_lim_cols = ['nimax', 'nimin']
for c in ni_lim_cols:
    if c not in df.columns:
        # If missing, create with defaults
        if c == 'nimax':
            df[c] = ni_max_default
        else:
            df[c] = 0

# Coerce numeric
df['nimax_num'] = pd.to_numeric(df['nimax'], errors='coerce').fillna(ni_max_default).astype(int)
df['nimin_num'] = pd.to_numeric(df['nimin'], errors='coerce').fillna(0).astype(int)

'''

# Global max and Big-M calcs based on ni limits
if specni:
    global_nimax = int(df['nimax'].max())
    sum_nimax_total = int(df['nimax'].sum())
else:
    global_nimax = ni_max_default
    sum_nimax_total = ni_max_default * len(df.index)

#Binary expension bits
#Ki_bits = math.ceil(math.log2(global_nimax + 1))

# Big-M upper bound for group-count-related constraints;
M_groups_val = min(most, sum_nimax_total)

##BUILDING OUR MODEL

#DEFINE SETS
model.g = Set(initialize=dfcp.index.unique().tolist()) #Groups and properties for GC Cp (Sahinidis et. al)
A_props=['Tci', 'Pci', 'Tbi',  'A0i', 'B0i', 'C0i', 'D0i', 'ai'] #To calculate Cp from Sahinidis
model.A=Set(initialize=A_props)

model.i = Set(initialize=df.index.unique().tolist()) #First order groups and properties for GC (Hukkerikar et. al)
X_props = ['Tb1i','Tm1i','δD1i','δP1i','δH1i','Vm1i']
model.X = Set(initialize=X_props) # The parameters calculated from MG method - might actually just be f(X) or not could put different constraints definition for some X 

model.B=Set(initialize=['acyclic','monocyclic']) #Molecule types; 'bicyclic' not included in this analysis

ni_lim=['nimax','nimin']
model.ni_lim = Set(initialize=ni_lim)

#Index sets for binary expansion
#model.Ki = RangeSet(1, Ki_bits)

#Assigning molecule types from Excel
TYPE_COL = "Type of molecule (like aromatic and so on) ?"
type_str_dict = df[TYPE_COL].to_dict()

model.Ga = Set(initialize=[i for i, t in type_str_dict.items() if t == 'aromatic'])
print("Ga groups:", ", ".join(sorted([str(i) for i in model.Ga])))
model.Gc = Set(initialize=[i for i, t in type_str_dict.items() if t == 'cyclic'])

#DEFINE PARAMETERS
if specni:
    nilim_data={(i, k): float(df.loc[i, k]) for i in df.index for k in ni_lim}
    model.nilim = Param(model.i, model.ni_lim, initialize=nilim_data, default=0)
    model.nilim.display()

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

# Big-M parameter for aromatic bonding
model.M_groups = Param(initialize=M_groups_val)

# Indicator Params (0/1) for membership and valency > 2
is_arom_dict = {i: 1 if i in list(model.Ga) else 0 for i in model.i}
model.is_arom = Param(model.i, initialize=is_arom_dict, default=0)

is_vgt2_arom_dict = {i: 1 if (i in list(model.Ga) and int(valency_dict.get(i, 0)) > 2) else 0 for i in model.i}
model.is_vgt2_arom = Param(model.i, initialize=is_vgt2_arom_dict, default=0)

#Mapping g groups to i groups - defined manually
ig_dict = {}
for i in dfSVSH.index:
    for g in dfSVSH.columns:
        try:
            ig_dict[(i, g)] = int(dfSVSH.loc[i, g]) if pd.notna(dfSVSH.loc[i, g]) else 0
        except:
            ig_dict[(i, g)] = 0
model.ig = Param(model.i, model.g, initialize=ig_dict, default=0)

#Constants
model.Tb0 = Param(initialize=244.7889) #Kelvins
model.Tm0 = Param(initialize=144.0977) #Kelvins
model.Vm0 = Param(initialize=0.0123) #m^3/kmol
model.R0 = Param(initialize=3.3) #MPa 1/2, based on experimental work 
model.sigmacD=Param(initialize=15.7) #MPa 1/2, CO2 solubility data
model.sigmacP=Param(initialize=5.2) #MPa 1/2
model.sigmacH=Param(initialize=5.8) #MPa 1/2
model.T_avg = Param(initialize=353) #in K, taken as average of absorption/desorption columns

#Target oriented Scaling parameters
model.t_rho = Param(initialize=17.4108)     # MEA target for rho scaling
model.t_Cp = Param(initialize=149.94)     # MEA target for Cp scaling
model.t_RED = Param(initialize=3.546) ## MEA target for RED scaling
model.d_rho = Param(initialize=5)     # allowed change for scaling
model.d_Cp = Param(initialize=30)     # allowed change for Cp scaling
model.d_RED = Param(initialize=3) ## allowed change for Cp scaling

# Objective function weights (mutable to allow updates in loop)
model.w_RED = Param(initialize=1.0, mutable=True)     # weight for RED
model.w_rho = Param(initialize=1.0, mutable=True)     # weight for density-NEGATIVE BECAUSE WE WANT TO MAXIMIZE IT
model.w_Cp  = Param(initialize=1.0, mutable=True)     # weight for heat capacity

#----DEFINE VARIABLES---## Need to justify bound selection
model.Tm=Var(within=NonNegativeReals, bounds=(0.01, 313), initialize=5) #in K
model.Tb=Var(within=NonNegativeReals, bounds=(393, 3000), initialize=500) #in K
model.rho=Var(within=NonNegativeReals, bounds=(0.0001, 50000)) #in mol.m^3
model.RED=Var(within=NonNegativeReals, bounds=(0.00001, None)) #in K
model.sol1=Var(within=Reals, bounds=(-16, 16)) #in 
model.sol2=Var(within=Reals, bounds=(-16, 16)) #in 
model.sol3=Var(within=Reals, bounds=(-16, 16)) #in 
model.molarvol=Var(within=NonNegativeReals, bounds=(0.000001, 100000)) #in 

#Set of binary variables (if we want to later define integar cuts to get rid of some solutions, also the variables are defined through linear expressions below)
#model.yi = Var(model.i, model.Ki, within=Binary)

#Binary molecule types : molecule type selection (acyclic vs monocyclic)
model.yb=Var(model.B, within=Binary) #binary molecule types defined later
model.m=Var(within=Reals, bounds=(-1,2))

#aromatic or cyclic mode 
model.ya =Var(within=Binary)  # aromatic mode Param (initialize=0) 
model.yc =Var(within=Binary)  # cyclic (non-aromatic) mode Param (initialize=0) 

#Integer counts (relaxed but integral via binary expansion)
model.ni = Var(model.i, within=NonNegativeIntegers, bounds=(0, global_nimax))
model.ng = Var(model.g, within=NonNegativeReals, bounds=(0, None), initialize=1.0)#could leave continuous potench

#Intermediate variables for calculating Cp, bounds copied from Sahinidis et. al where bounding doesn't result in infeasible
#Could eliminate by defining Cp but then can't bound
model.T_b = Var(within=NonNegativeReals, initialize=350, bounds=(50, 1000)) #K
model.T_c = Var(within=NonNegativeReals, initialize=500, bounds=(100, 2000)) #K
model.P_c = Var(within=PositiveReals, initialize=5.0, bounds=(2, 200)) #bar
model.T_avgr = Var(within=PositiveReals, initialize=0.7, bounds=(0.001, 1)) #K/K
model.alpha = Var(within=Reals, initialize=0.1, bounds=(-5.4, 6100)) # N/A
model.beta = Var(within=Reals, initialize=-1.0, bounds=(-160000, 0)) # N/A
model.W = Var(within=Reals, initialize=0.3, bounds=(-1, 1.3)) # N/A
model.Cp0 = Var(within=NonNegativeReals, initialize=1) #cal/mol*K (converted so Cp is in J/mol*K)

#Properties results
model.Cp = Var(within=NonNegativeReals, bounds=(10e-6, 1000)) #J/mol*K
model.fX = Var(model.X, within=Reals)  #Tm, Tb, etc., Might be an issue that we can't/haven't initialized them independently

# Binary variable indicating whether attachment out of the aromatic ring is permitted
model.z_attach_ok = Var(within=Binary)

# Expressions for aromatic counts (should be in constraints but since not directly the constraint ok here)
model.count_arom_vgt2 = Expression(expr=sum(model.ni[i] * model.is_vgt2_arom[i] for i in model.i))
model.count_non_aromatic = Expression(expr=sum(model.ni[i] * (1 - model.is_arom[i]) for i in model.i))

##DEFINE CONSTRAINTS
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
    return model.T_c * (0.584 + 0.965 * sum(model.ng[g] * model.cg[g, 'Tci'] for g in model.g) - (sum(model.ng[g] * model.cg[g, 'Tci'] for g in model.g))**2 ) == model.T_b
model.Tc_calc = Constraint(rule=Tc_calc_rule)

# Critical pressure
def Pc_calc_rule(model):
    return model.P_c == 1 / (0.113 + 0.0032 * sum(model.ng[g] * model.cg[g, 'ai'] for g in model.g) - sum(model.ng[g] * model.cg[g, 'Pci'] for g in model.g))**(2)
model.Pc_calc = Constraint(rule=Pc_calc_rule)

# Reduced average temperature (reformulated)
def Tavgr_calc_rule(model):
    return model.T_avgr * model.T_c == model.T_avg
model.Tavgr_calc = Constraint(rule=Tavgr_calc_rule)

# Acentric factor
# Split alpha and beta into rules because somehow it changes the calculation??
def alpha_calc_rule(model):
    return model.alpha == (-5.97214 - log(model.P_c / 1.013) + ((6.09648 * model.T_c) / model.T_b) + 
             1.28862 * log(model.T_b / model.T_c) - 0.169347 * (model.T_b / model.T_c)**6)
model.alpha_calc = Constraint(rule=alpha_calc_rule)

def beta_calc_rule(model):
    return model.beta == (15.2518 - ((15.6875 * model.T_c) / model.T_b) - 13.4721 * log(model.T_b / model.T_c) + 
            0.43577 * (model.T_b / model.T_c)**6)
model.beta_calc = Constraint(rule=beta_calc_rule)

def W_calc_rule(model):
    #model.alpha = (-5.97214 - log((model.P_c / 1.013) + (6.09648 * model.T_c) / model.T_b) + 
    #         1.28862 * log(model.T_b / model.T_c) - 0.167347 * (model.T_b / model.T_c)**6)
    #model.beta = (15.2518 - ((15.6875 * model.T_c) / model.T_b) - 13.4721 * log(model.T_b / model.T_c) + 
    #        0.43577 * (model.T_b / model.T_c)**6)
    return model.W * model.beta == model.alpha
model.W_calc = Constraint(rule=W_calc_rule)

# Ideal liquid heat capacity
def Cp0_calc_rule(model):
    return model.Cp0 == (sum(model.ng[g] * model.cg[g, 'A0i'] for g in model.g) - 37.93 +
                         (sum(model.ng[g] * model.cg[g, 'B0i'] for g in model.g) + 0.21) * model.T_avg +
                         (sum(model.ng[g] * model.cg[g, 'C0i'] for g in model.g) - 3.91e-4) * model.T_avg**2 +
                         (sum(model.ng[g] * model.cg[g, 'D0i'] for g in model.g) + 2.06e-7) * model.T_avg**3)
model.Cp0_calc = Constraint(rule=Cp0_calc_rule)

# Final Cp liquid calculation (Rowlinson)
def Cp_rule(model):#either we define specific exceptions to take into account the various f(X) or idk
    return model.Cp == 4.184 * ((1/4.1868)*(model.Cp0 + 8.314*(1.45 + (0.45/(1-model.T_avgr)) + 0.25*model.W * (17.11 + 25.2*(((1 - model.T_avgr)**(1/3))/model.T_avgr) + (1.742 / (1-model.T_avgr))))))
model.Cp_constraint = Constraint(rule=Cp_rule)

# (4) Link Hurekkikar group contributions to actual properties
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

# (5) Definition of density from GC calculated parameters
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

#(9) Expressing ni as an integer variables with linear expressions - NOT USED ATM
#def ni_rule(model,i):
#    return model.ni[i] == sum(2**(k-1)*model.yi[i,k] for k in model.Ki)
#model.ni_rule = Constraint(model.i,rule=ni_rule)

#(11) Octet Rule, i compounds
def vi_rule(model):
    return 2*model.m==sum((2-model.vi[i])*model.ni[i] for i in model.i)
model.vi_rule = Constraint(rule=vi_rule)

#(12) One cyclic mode
# Exactly one mode if monocyclic, none if acyclic
def monocyclic_mode_selection_rule(model):
    return model.ya + model.yc == model.yb['monocyclic']
model.monocyclic_mode_selection = Constraint(rule=monocyclic_mode_selection_rule)

#(13) Aromatic & cyclic constraints
# Exactly 6 aromatic groups if aromatic mode; 0 aromatic groups otherwise -- SHOULD BE 5 OR 6
def aromatic_ring_rule(model):
    return sum(model.ni[i] for i in model.Ga) == 6 * model.ya
model.aromatic_ring = Constraint(rule=aromatic_ring_rule)

# Forbid cyclic groups unless cyclic mode is chosen (Big-M)- Check If getting rid of this still work
def cyclic_only_if_cyclic_mode_rule(model):
    return sum(model.ni[i] for i in model.Gc) <= model.M_groups * model.yc
#model.cyclic_only_if_cyclic_mode = Constraint(rule=cyclic_only_if_cyclic_mode_rule)

def cyclic_upper_only_if_cyclic_mode_rule(model):
    return sum(model.ni[i] for i in model.Gc) <= 8 * model.yc
model.cyclic_only_if_cyclic_mode = Constraint(rule=cyclic_upper_only_if_cyclic_mode_rule)

def cyclic_lower_only_if_cyclic_mode_rule(model):
    return sum(model.ni[i] for i in model.Gc) >= 5 * model.yc
model.cyclic_only_if_cyclic_mode = Constraint(rule=cyclic_lower_only_if_cyclic_mode_rule)

#(14) Aromatic group with a valency of >2 only can bind to non-aromatic groups
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
    return sum(model.ni[k] for k in model.i)>=least
model.zer_rule = Constraint(rule=zer)    

#(17) Maximum amount of groups
def maxgroups(model):
    return sum(model.ni[k] for k in model.i)<=most
model.max_rule = Constraint(rule=maxgroups)    

#(18) Minimum cyclic groups if cyclic non-aromatic mode : add a max cyclic ?
def min_cyclic_groups_rule(model):
    return sum(model.ni[i] for i in model.Gc) >= 5 * model.yc
model.min_cyclic_groups = Constraint(rule=min_cyclic_groups_rule)


#(19) Min and max ni constraints
if specni:
        print("adding nimin/max constraints")
        def minmaxni_rule(model, i, k):
            if k == 'nimax':
                return model.ni[i] <= model.nilim[i, 'nimax']
            elif k == 'nimin':
                return model.ni[i] >= model.nilim[i, 'nimin']
        model.nibounds = Constraint(model.i, model.ni_lim, rule=minmaxni_rule)

##OBJECTIVE FUNCTION
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

weights = [[1, 1, 1], [0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 0]]

summary_rows = []
ni_rows = []
ng_rows = []

# Prepare Excel writer to same folder as input
out_dir = os.path.dirname(inputGCs) if os.path.dirname(inputGCs) else os.getcwd()
out_path = os.path.join(out_dir, "APO_results.xlsx")

#helper function for feasibility data
def compute_feasibility_stats(model, tol=1e-6):
    max_constr_viol = 0.0
    violated_constrs = 0
    total_constrs = 0

    for c in model.component_objects(Constraint, active=True):
        for _, ci in c.items():
            total_constrs += 1
            body = value(ci.body, exception=False)
            if body is None:
                # Cannot evaluate now; skip
                continue
            lb = value(ci.lower, exception=False)
            ub = value(ci.upper, exception=False)
            viol = 0.0
            if lb is not None:
                viol = max(viol, max(0.0, lb - body))
            if ub is not None:
                viol = max(viol, max(0.0, body - ub))
            if viol > tol:
                violated_constrs += 1
            if viol > max_constr_viol:
                max_constr_viol = viol

    max_var_viol = 0.0
    for v in model.component_data_objects(Var, active=True):
        val = value(v, exception=False)
        if val is None:
            continue
        lb = v.lb
        ub = v.ub
        viol = 0.0
        if lb is not None:
            viol = max(viol, max(0.0, lb - val))
        if ub is not None:
            viol = max(viol, max(0.0, val - ub))
        if viol > max_var_viol:
            max_var_viol = viol

    return dict(
        max_constr_viol=max_constr_viol,
        max_var_viol=max_var_viol,
        violated_constrs=violated_constrs,
        total_constrs=total_constrs
    )

def safe_val(expr):
    return value(expr, exception=False)

#SOLVER OPTIONS
solver = SolverFactory('gams')
solver.options['solver'] = 'DICOPT'

for a, current_weights in enumerate(weights):
    # Update mutable Params
    model.w_RED.set_value(current_weights[0])
    model.w_rho.set_value(current_weights[1])
    model.w_Cp.set_value(current_weights[2])

    # Solve
    results = solver.solve(model, tee=False)

    # Solver/termination info
    solver_status = results.solver.status if hasattr(results, 'solver') else None
    term_cond = results.solver.termination_condition if hasattr(results, 'solver') else None

    # Compute feasibility metrics
    feas_stats = compute_feasibility_stats(model, tol=1e-6)
    max_constr_viol = feas_stats['max_constr_viol']
    max_var_viol = feas_stats['max_var_viol']

    # Define feasibility flag:
    # True if small violations AND termination suggests feasibility/optimality
    tc = term_cond
    tc_ok = tc in (TerminationCondition.optimal,
                   TerminationCondition.locallyOptimal,
                   TerminationCondition.feasible)
    small_viol = (max_constr_viol <= 1e-6) and (max_var_viol <= 1e-6)
    is_feasible = bool(tc_ok and small_viol)

    # Extract values safely
    objv = safe_val(model.obj)
    cpv = safe_val(model.Cp)
    rhov = safe_val(model.rho)
    redv = safe_val(model.RED)
    tmv = safe_val(model.Tm)
    tbv = safe_val(model.Tb)
    vv = safe_val(model.molarvol)

    y_acyclic = safe_val(model.yb['acyclic']) if 'acyclic' in model.B else 0
    y_monocyclic = safe_val(model.yb['monocyclic']) if 'monocyclic' in model.B else 0
    m_type = 'monocyclic' if y_monocyclic and y_monocyclic > 0.5 else 'acyclic'
    cyc_mode = 'aromatic' if (safe_val(model.ya) and safe_val(model.ya) > 0.5) else (
        'non-aromatic cyclic' if (safe_val(model.yc) and safe_val(model.yc) > 0.5) else 'none'
    )

    summary_rows.append({
        'run': a,
        'w_RED': current_weights[0],
        'w_rho': current_weights[1],
        'w_Cp': current_weights[2],
        'objective': objv,
        'Cp_J_per_molK': cpv,
        'rho_mol_per_m3': rhov,
        'RED': redv,
        'Tm_K': tmv,
        'Tb_K': tbv,
        'molarvol_m3_per_kmol': vv,
        'm_type': m_type,
        'cyclic_mode': cyc_mode,
        'solver_status': str(solver_status),
        'termination_condition': str(term_cond),
        'is_feasible': is_feasible,
        'max_constr_viol': max_constr_viol,
        'max_var_viol': max_var_viol
    })

    for ii in model.i:
        niv = safe_val(model.ni[ii])
        if niv is not None and niv > 1e-6:
            ni_rows.append({'run': a, 'group_i': str(ii), 'ni': niv})

    for g in model.g:
        ngv = safe_val(model.ng[g])
        if ngv is not None and ngv > 1e-6:
            ng_rows.append({'run': a, 'group_g': str(g), 'ng': ngv})

# Write Excel report
with pd.ExcelWriter(out_path, engine='xlsxwriter') as writer:
    pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)
    pd.DataFrame(ni_rows).to_excel(writer, sheet_name='ni_counts', index=False)
    pd.DataFrame(ng_rows).to_excel(writer, sheet_name='ng_counts', index=False)

print(f"\nResults written to: {out_path}")

# ---- DISPLAY RESULTS (last run) ----
print("\n" + "="*70)
print("OPTIMAL SOLUTION FOR CARBON CAPTURE SOLVENT DESIGN")
print("="*70)
print(f"weights used (last run): {weights[-1]}")

print("\nGroup Counts (ni) - Hukkerikar Groups:")
print("-" * 70)
for ii in model.i:
    if model.ni[ii].value and model.ni[ii].value > 0.01:
        print(f"  {ii:25s}: {model.ni[ii].value:6.3f}")

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
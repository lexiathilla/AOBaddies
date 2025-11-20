from pyomo.environ import *
import math
import matplotlib.pyplot as plt
model = ConcreteModel('APOProject')

#max number of groups per elements
ni_max=8
ng_max=10

#Define SETS
model.g=Set(initialize=['CH3', 'CH2','CH','C','CHdouble','Cdouble','Cdouble','NH2','NH','N','Ndouble','OH','O','Odouble'])#groups for GC Cp
model.i=Set(initialize=['CH3','CH2','CH','C','CH2doubleCH','CHdoubleCH','CH2doubleC'])#First order groups for GC - need to choose from large variety - litterature backed ?
model.A=Set(initialize=['Cp'])
model.X=Set(initialize=['Tm','Tb','molarvol','sol1','sol2','sol3']) #The parameters calculated from MG method - might actually just be f(X) or not could put different constraints definition for some X 
model.B=Set(initialize=['acyclic','monocyclic','bicyclic'])

#Ki={i: math.ceil(math.log2(ni)) for i in model.i}  could create ngs and nis for each g and i btw
#Kg={g: math.ceil(math.log2(ng)) for g in model.g}
Ki = math.ceil(math.log2(ni_max))
Kg = math.ceil(math.log2(ng_max))
# Create index sets
model.Ki = Set(1, Ki)
model.Kg = Set(1, Kg)

#Define Parameters GC, ideally for each grp considered and each property (i,X)
i_data={('CH3', 'Tm'): 1, ('CH2', 'Tm'): 2, ('CH', 'Tm'): 2, ('C', 'Tm'): 1, ('CH2doubleCH', 'Tm'): 1,('CHdoubleCH', 'Tm'): 1,('CH2doubleC', 'Tm'): 1
        ,('CH3', 'Tb'): 1, ('CH2', 'Tb'): 2, ('CH', 'Tb'): 2, ('C', 'Tb'): 1, ('CH2doubleCH', 'Tb'): 1,('CHdoubleCH', 'Tb'): 1,('CH2doubleC', 'Tb'): 1
        ,('CH3', 'molarvol'): 1, ('CH2', 'molarvol'): 2, ('CH', 'molarvol'): 2, ('C', 'molarvol'): 1,('CH2doubleCH', 'molarvol'): 1,('CHdoubleCH', 'molarvol'): 1,('CH2doubleC', 'molarvol'): 1
        ,('CH3', 'sol1'): 1, ('CH2', 'sol1'): 2, ('CH', 'sol1'): 2, ('C', 'sol1'): 1, ('CH2doubleCH', 'sol1'): 1,('CHdoubleCH', 'sol1'): 1,('CH2doubleC', 'sol1'): 1
        ,('CH3', 'sol2'): 1, ('CH2', 'sol2'): 2, ('CH', 'sol2'): 2, ('C', 'sol2'): 1, ('CH2doubleCH', 'sol2'): 1,('CHdoubleCH', 'sol2'): 1,('CH2doubleC', 'sol2'): 1
        ,('CH3', 'sol3'): 1, ('CH2', 'sol3'): 2, ('CH', 'sol3'): 2, ('C', 'sol3'): 1,('CH2doubleCH', 'sol3'): 1,('CHdoubleCH', 'sol3'): 1,('CH2doubleC', 'sol3'): 1}
model.ci=Param(model.i, model.X, initialize=i_data, default=0)

#for Cp
g_data={('CH3', 'Cp'): 1, ('CH2', 'Cp'): 2, ('CH', 'Cp'): 2, ('C', 'Cp'): 1, ('CHdouble', 'Cp'): 1,('Cdouble', 'Cp'): 1,('NH2', 'Cp'): 1,('NH', 'Cp'): 1,('N', 'Cp'): 1,('Ndouble', 'Cp'): 1,('OH', 'Cp'): 1,('O', 'Cp'): 1,('Odouble', 'Cp'): 1 }
model.cg=Param(model.g, model.A, initialize=g_data, default=0)

#valency data (defined for i for now)
v_i={'CH3':1, 'CH2':2, 'CH':3, 'C':4, 'CH2doubleCH':1, 'CHdoubleCH':3, 'CH2doubleC':2}
model.vi=Param(model.i, initialize=v_i, default=0)
v_g={}
model.vg=Param(model.i, initialize=v_g, default=0)
#group types : 1=general, 2=unsaturated, 3=aromatic; probably wrong way to do it 
type_data={'CH3':1, 'CH2':1, 'CH':1, 'C':1, 'CH2doubleCH':2, 'CHdoubleCH':2, 'CH2doubleC':2}
model.grouptype = Param(model.i, initialize=type_data, default=0)

#Number of groups in other groups (link between Cp and Hurekkikar et al) #TO CORRECT, when relation not mentionned we get zero
groups_data={('CH3', 'CH3'): 1,
             ('CH2', 'CH2'): 1, 
             ('CH', 'CH'): 1, 
             ('C', 'C'): 1, 
             ('CH2doubleCH', 'CHdouble'): 1,('CH2doubleCH', 'CH2'): 1,
             ('CHdoubleCH', 'CHdouble'): 2,
             ('CH2doubleC', 'CH2'): 1,('CH2doubleC', 'Cdouble'): 1 }
model.ig=Param(model.i, model.g, initialize=groups_data, default=0)

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
#Binary molecule types
model.yb=Var(model.B, within=Binary) #binary molecule types defined later
model.m=Var(within=Reals)
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
    return model.fX[X]==(sum(model.ni[i] * model.ci[i]) for i in model.i)
model.XGC=Constraint(model.X,rule=XGC_rule)

# (2) Number of g-groups derived from i-groups
def numgini_rule(model,g):#either we define specific exceptions to take into account the various f(X) or idk
    return model.ng[g]==sum(model.ni[i] * model.ig[i,g] for i in model.i)
model.numgini = Constraint(model.g, rule=numgini_rule)

# (3) Cp calculations from g-groups
def CpGC_rule(model):#either we define specific exceptions to take into account the various f(X) or idk
    return model.Cp == sum(model.ng[g] * model.cg[g, 'Cp'] for g in model.g)
model.CpGC = Constraint(rule=CpGC_rule)

# (4) def the various functions in Hurekkikar group contributions
def link_properties_rule(model, X):
    if X == 'Tm':
        return model.Tm == exp(model.fX[X]/model.Tm0)
    elif X == 'Tb':
        return model.Tb == exp(model.fX[X]/model.Tb0)
    elif X == 'molarvol':
        return model.molarvol == model.fX[X]+model.Vm0
    elif X in ['sol1', 'sol2', 'sol3']:
        return getattr(model,X) == model.fX[X]
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
    return model.m == model.yb['acyclic'] - model.yb['bicyclic']
model.m_rule = Constraint(rule=m_rule)

#(9) Replacing integer variables with linear expressions : i 
def ni_rule(model,i):
    return model.ni[i] == sum(2**(k-1)*model.yi[i,k] for k in model.Ki)
model.ni_rule = Constraint(rule=ni_rule)

#(10) Replacing integer variables with linear expressions : g
def ng_rule(model,g):
    return model.ng[g] == sum(2**(k-1)*model.yg[g,k] for k in model.Kg)
model.ng_rule = Constraint(rule=ng_rule)

#(11) Octet Rule, i compounds
def vi_rule(model):
    return 2*model.m==sum((2-model.vi[i])*model.ni[i] for i in model.i)
model.vi_rule = Constraint(rule=vi_rule)

#(12) Octet Rule, g compounds
def vg_rule(model):
    return 2*model.m==sum((2-model.vg[g])*model.ng[g] for g in model.g)
model.vg_rule = Constraint(rule=vg_rule)

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
solver=SolverFactory('gams:conopt')
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


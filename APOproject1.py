from pyomo.environ import *
import math
import matplotlib.pyplot as plt
model = ConcreteModel('APOProject')

"""
yadayada"""

#Define SETS
#groups for GC Cp -- from Sahinidis method, which uses Joback & Reid groups
    #only need to include what's needed to build chosen i groups
model.g=Set(initialize=['CH3', 'CH2','CH','C','CH2d','CHd','Cd','Cdd','CHt','Ct',               #Non-ring increments
                        'CH2r', 'CHr', 'Cr', 'CHdr', 'Cdr',                                     #Ring increments
                        'F', 'Cl', 'Br', 'I',                                                   #Halogen increments 
                        'OHalc', 'OHphe', 'O', 'Or', 'CO', 'COr', 'OCH', 'COOH', 'COO', 'Od',   #Oxygen increments (ring and non-ring)
                        'NH2','NH', 'NHr', 'N', 'Nd', 'Ndr', 'NHd', 'CN', 'NO2',                #Nitrogen increments
                        'SH', 'S', 'Sr'])                                                       #Sulfur increments
model.A=Set(initialize=['Cp']) #Calculated from Sahinidis

#First order groups for GC from Hukkerikar - need to finalise selection
model.i=Set(initialize=['CH3','CH2','CH','C','CH2doubleCH','CHdoubleCH','CH2doubleC', 'CHdoubleC', 'CdoubleC', 'CH2doubleCdoubleCH'])#...
model.X=Set(initialize=['Tm','Tb','molarvol','sol1','sol2','sol3']) #The parameters calculated from MG method - might actually just be f(X) or not could put different constraints definition for some X 
model.B=Set(initialize=['acyclic','monocyclic','bicyclic'])

#Define Parameters GC, ideally for each grp considered and each property (i,X)
i_data={('CH3', 'Tm'): 1, ('CH2', 'Tm'): 2, ('CH', 'Tm'): 2, ('C', 'Tm'): 1, ('CH2doubleCH', 'Tm'): 1,('CHdoubleCH', 'Tm'): 1,('CH2doubleC', 'Tm'): 1
        ,('CH3', 'Tb'): 1, ('CH2', 'Tb'): 2, ('CH', 'Tb'): 2, ('C', 'Tb'): 1, ('CH2doubleCH', 'Tb'): 1,('CHdoubleCH', 'Tb'): 1,('CH2doubleC', 'Tb'): 1
        ,('CH3', 'molarvol'): 1, ('CH2', 'molarvol'): 2, ('CH', 'molarvol'): 2, ('C', 'molarvol'): 1,('CH2doubleCH', 'molarvol'): 1,('CHdoubleCH', 'molarvol'): 1,('CH2doubleC', 'molarvol'): 1
        ,('CH3', 'sol1'): 1, ('CH2', 'sol1'): 2, ('CH', 'sol1'): 2, ('C', 'sol1'): 1, ('CH2doubleCH', 'sol1'): 1,('CHdoubleCH', 'sol1'): 1,('CH2doubleC', 'sol1'): 1
        ,('CH3', 'sol2'): 1, ('CH2', 'sol2'): 2, ('CH', 'sol2'): 2, ('C', 'sol2'): 1, ('CH2doubleCH', 'sol2'): 1,('CHdoubleCH', 'sol2'): 1,('CH2doubleC', 'sol2'): 1
        ,('CH3', 'sol3'): 1, ('CH2', 'sol3'): 2, ('CH', 'sol3'): 2, ('C', 'sol3'): 1,('CH2doubleCH', 'sol3'): 1,('CHdoubleCH', 'sol3'): 1,('CH2doubleC', 'sol3'): 1}
model.ci=Param(model.i, model.X, initialize=i_data, default=0)

#for Cp
#Have only defined for non-ring increments to test - require Pci, Tci, Tbi, Cp0i(a,b,c,d) for each group
#Unsure if Pci, Tci, Tbi compatible with Hukkerikar values? May not need to define twice 
g_data={('CH3', 'Tci'):0.0141, ('CH2', 'Tci'):0.0189, ('CH', 'Tci'):0.0164, ('C', 'Tci'):0.0067, ('CH2d', 'Tci'):0.0113, ('CHd', 'Tci'):0.0129, ('Cd', 'Tci'):0.0117, ('Cdd', 'Tci'):0.0026, ('CHt', 'Tci'):0.0027, ('Ct', 'Tci'):0.0020,
        ('CH3', 'Tbi'):23.58, ('CH2', 'Tbi'):22.88, ('CH', 'Tbi'):21.74, ('C', 'Tbi'):18.25, ('CH2d', 'Tbi'):18.18, ('CHd', 'Tbi'):24.96, ('Cd', 'Tbi'):24.14, ('Cdd', 'Tbi'):26.15, ('CHt', 'Tbi'):9.20, ('Ct', 'Tbi'):27.38,
        ('CH3', 'Pci'):-0.0012, ('CH2', 'Pci'):0, ('CH', 'Pci'):0.0020 ('C', 'Pci'):0.0043, ('CH2d', 'Pci'):-0.0028, ('CHd', 'Pci'):-0.0006, ('Cd', 'Pci'):0.0011, ('Cdd', 'Pci'):0.0028, ('CHt', 'Pci'):-0.0008, ('Ct', 'Pci'):0.0016,
        #NOT FILLED IN
        ('CH3', 'A0i'):0.0141, ('CH2', 'A0i'):0.0189, ('CH', 'A0i'):0.0164, ('C', 'A0i'):0.0067, ('CH2d', 'A0i'):0.0113, ('CHd', 'A0i'):0.0129, ('Cd', 'A0i'):0.0117, ('Cdd', 'A0i'):0.0026, ('CHt', 'A0i'):0.0027, ('Ct', 'A0i'):0.0020,
        ('CH3', 'B0i'):0.0141, ('CH2', 'B0i'):0.0189, ('CH', 'B0i'):0.0164, ('C', 'B0i'):0.0067, ('CH2d', 'B0i'):0.0113, ('CHd', 'B0i'):0.0129, ('Cd', 'B0i'):0.0117, ('Cdd', 'B0i'):0.0026, ('CHt', 'B0i'):0.0027, ('Ct', 'B0i'):0.0020,
        ('CH3', 'C0i'):0.0141, ('CH2', 'C0i'):0.0189, ('CH', 'C0i'):0.0164, ('C', 'C0i'):0.0067, ('CH2d', 'C0i'):0.0113, ('CHd', 'C0i'):0.0129, ('Cd', 'C0i'):0.0117, ('Cdd', 'C0i'):0.0026, ('CHt', 'C0i'):0.0027, ('Ct', 'C0i'):0.0020,
        ('CH3', 'D0i'):0.0141, ('CH2', 'D0i'):0.0189, ('CH', 'D0i'):0.0164, ('C', 'D0i'):0.0067, ('CH2d', 'D0i'):0.0113, ('CHd', 'D0i'):0.0129, ('Cd', 'D0i'):0.0117, ('Cdd', 'D0i'):0.0026, ('CHt', 'D0i'):0.0027, ('Ct', 'D0i'):0.0020,
        }
model.cg=Param(model.g, model.A, initialize=g_data, default=0)

#valency data (defined for i for now)
v_i={'CH3':1, 'CH2':2, 'CH':3, 'C':4, 'CH2doubleCH':1, 'CHdoubleCH':3, 'CH2doubleC':2}
model.vi=Param(model.i, initialize=v_i, default=0)

#group types : 1=general, 2=unsaturated, 3=aromatic; probably wrong way to do it 
type_data={'CH3':1, 'CH2':1, 'CH':1, 'C':1, 'CH2doubleCH':2, 'CHdoubleCH':2, 'CH2doubleC':2}
model.grouptype = Param(model.i, initialize=type_data, default=0)

#Number of groups in other groups (link between Cp and Hurekkikar et al)
#Have defined for 1-10 Hukkerikar groups TO TEST 
# #TO CORRECT, when relation not mentioned we get zero
groups_data={('CH3', 'CH3'): 1,
             ('CH2', 'CH2'): 1, 
             ('CH', 'CH'): 1, 
             ('C', 'C'): 1, 
             ('CH2doubleCH', 'CHd'): 1,('CH2doubleCH', 'CH2'): 1,
             ('CHdoubleCH', 'CHd'): 2,
             ('CH2doubleC', 'CH2'): 1,('CH2doubleC', 'Cd'): 1 
             ('CdoubleC', 'Cd'): 2
             ('CHdoubleCdoubleCH', 'CH2d'):1, ('CHdoubleCdoubleCH', 'Cd'):1,('CHdoubleCdoubleCH', 'CHd'):1}
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

#Binary molecule types
model.yb=Var(model.B, within=Binary)
model.m=Var(within=Reals)
#Integer count (to be replaced with binary variables)
model.ni=Var(model.i, within=NonNegativeInteger)#!!! watch w the integer cuts thing might have to def ni=sumk=0 to k 2^kyik (continuous associated with a set of binary variables)
model.ng = Var(model.g, within=NonNegativeInteger)

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
        return model.Tm == log(model.fX[X])*model.Tm0
    elif X == 'Tb':
        return model.Tb == log(model.fX[X])*model.Tb0
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

#(9) Indicator variables to numbers of group
#Octet rule
def rho_def_rule(model):
    return model.rho == 1 / model.molarvol
model.rho_def = Constraint(rule=rho_def_rule)
#

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


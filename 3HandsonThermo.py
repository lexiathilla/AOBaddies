from pyomo.environ import *
import math
model = ConcreteModel('Thermo2')

#minimizing free energy through det eqm compo of mixture 0.5N2H4 (A)+ 0.502 (B)
#Defining Parameters
model.P=750#psi /14.696 to get atm
model.T=3500#K
##Defining variables
#model.A = Var(bounds=(0,1))
#model.B = Var(bounds=(0,1))
#OR alternatively use indexed comp for tractability 
#components = ['A', 'B']
#model.x = Var(components, bounds=(0, 1))
#AND NO ! eheh take into account all elements in a set 

#Define SETS
model.i = Set(initialize=['N','H','O']) #Elements (3)
model.c = Set(initialize=['N','N2','NO','NH','H','H2','H2O','O','O2','OH']) #Compounds (10)

#Define Parameters
#initial number of each atoms 0.5N2H4 (A)+ 0.502 (B)
b_data={'N':1,'H':2,'O':1}
model.b=Param(model.i, initialize=b_data)

#ATOM COMPOSITION TABLE a(i,c) but not in the same order
a_data = {
    ('H', 'H'): 1, ('H', 'H2'): 2, ('H', 'H2O'): 2, ('H', 'NH'): 1, ('H', 'OH'): 1,
    ('N', 'N'): 1, ('N', 'N2'): 2, ('N', 'NH'): 1, ('N', 'NO'): 1,
    ('O', 'H2O'): 1, ('O', 'NO'): 1, ('O', 'O'): 1, ('O', 'O2'): 2, ('O', 'OH'): 1}

# Flip (compound, element) -> (element, compound)
a_fixed = {(elem, comp): val for (comp, elem), val in a_data.items()}
model.a=Param(model.i,model.c, initialize=a_fixed, default=0)#if keys spec., no need to have table in order

#GIBBS FREE ENERGY STUFF for individual compounds
g_data={'H':-10.021,'H2':-21.096,'H2O':-37.986,'N':-9.846,'N2':-28.653,'NH':-18.918,'NO':-28.032,'O':-14.640, 'O2':-30.594,'OH':-26.11}

model.g=Param(model.c,initialize=g_data)

#pressure correction
gibbs_p_data={c: g_data[c]+math.log(model.P*0.07031) for c in g_data.keys()}
model.gp=Param(model.c, initialize=gibbs_p_data)

#DEFINE VARIABLES
model.xc = Var(model.c, within=NonNegativeReals, bounds=(0.001, None)) 
model.xt = Var(within=NonNegativeReals, bounds=(0.01, None))

# --- Constraints ---
def total_moles_rule(model):
    return model.xt == sum(model.xc[c] for c in model.c)
model.total_moles = Constraint(rule=total_moles_rule)

def balance_rule(model, i):#because done for all i 
    return sum(model.a[i,c]*model.xc[c] for c in model.c)==model.b[i]
model.balance = Constraint(model.i, rule=balance_rule)

#Objective : min tot gibbs free energy

def gibbs_free_energy(model):
    return sum(model.xc[c]*(model.gp[c])+log(model.xc[c]/(model.xt)) for c in model.c) #!!!don't use math.log on pyomo variables that python cannot convert to float
model.obj=Objective(rule=gibbs_free_energy, sense=minimize)

#Solve
solver=SolverFactory('gams:conopt')
results=solver.solve(model, tee=True)

#Display results
# --- Display results ---
print("\n--- Optimal Solution ---")
for c in model.c:
    print(f"xc[{c}] = {value(model.xc[c]):.6f}")
print(f"\nx_tot = {value(model.xt):.6f}")
print(f"G = {value(model.obj):.6f}")

# --- Optional: plot ---

import matplotlib.pyplot as plt
x_vals = [value(model.xc[c]) for c in model.c]
mole_frac = [x / value(model.xt) for x in x_vals]

plt.figure(figsize=(8,5))
plt.bar(model.c, mole_frac)
plt.title("Equilibrium Mole Fraction per Compound")
plt.ylabel("Mole Fraction")
plt.xticks(rotation=45)
plt.grid(axis='y', linestyle='--', alpha=0.6)
plt.show()

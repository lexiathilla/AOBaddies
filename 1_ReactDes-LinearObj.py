
from pyomo.environ import *

# Create a concrete model
model = ConcreteModel('ReactDes')

# Sets
model.I = Set(initialize=['Ain', 'A', 'B'])
model.J = Set(initialize=[1, 2])

# Parameters
model.k = Param(model.J, initialize={1: 5, 2: 2}) # Reaction rate
model.V = Param(initialize=1)  # Reactor Volume

# Variables
model.C = Var(model.I, domain=NonNegativeReals)
model.F = Var(domain=NonNegativeReals)

# Objective
def obj_rule(m):
    return m.C['B']
model.obj = Objective(expr=model.C['B'], sense=maximize)

# Equations / Constraints
def eq1_rule(m):
    return m.F * (m.C['Ain'] - m.C['A']) - m.k[1] * m.C['A'] * m.V == 0
model.eq1 = Constraint(rule=eq1_rule)

def eq2_rule(m):
    return -m.F * m.C['B'] + (m.k[1] * m.C['A'] - m.k[2] * m.C['B']) * m.V == 0
model.eq2 = Constraint(rule=eq2_rule)

# Bounds
model.C['Ain'].setlb(0)
model.C['Ain'].setub(1)
model.C['A'].setub(0.1)
model.F.setlb(0)
model.F.setub(20)

# Initial values
model.C['A'].value = 0.01
model.C['Ain'].value = 0.8
model.C['B'].value = 0.5
model.F.value = 10

model.C['Ain'].value=0.1

# Solve
solver = SolverFactory('gams:baron')
solver.solve(model, tee=True)

# Display all results
model.display()



from pyomo.environ import *

# --- Create model ---
model = ConcreteModel('Gardener')

# --- Parameters ---
model.F  = Param(initialize=16)                    # fencing available

# --- Variables ---
model.L = Var(within=NonNegativeReals)
model.W = Var(within=NonNegativeReals)
model.A = Var(within=NonNegativeReals)

# --- Constraints ---
#! there's a certain number of supported nonlinear expressions
#  write constraints
model.area = Constraint(expr=model.A==model.L*model.W)
model.perimeter = Constraint(expr= model.F==2*(model.L+model.W))

# --- Objective (maximize area) ---
#  write objective
model.obj=Objective(expr=model.A, sense=maximize)
# --- Initialize ---
W0 = 1                 # initial width
L0 = 7                    # initial length
A0 = W0 * L0  # initial area
#  write initialization
model.L.value = L0
model.W.value =W0
model.A.value =A0
# --- Solve ---
solver = SolverFactory('gams:conopt')   # nonlinear solver
results = solver.solve(model, tee=True)

# --- Display results ---
print(f"Optimal Length (L): {value(model.L):.3f}")
print(f"Optimal Width  (W): {value(model.W):.3f}")
print(f"Optimal Area   (A): {value(model.A):.3f}")

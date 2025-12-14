#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pyomo.environ import (
    ConcreteModel, Set, Param, Var, Binary,
    NonNegativeReals, NonNegativeIntegers, PositiveReals, Reals,
    Expression, Constraint, Objective, minimize, log, exp, sqrt, value, SolverFactory
)
from pyomo.opt import TerminationCondition
import pandas as pd
import math
import json

# ------------------------------------------------------------
# 1) Model builder: wrap your APOprojectint.py content here
# ------------------------------------------------------------
def build_model(xlsx_path: str) -> ConcreteModel:
    model = ConcreteModel("APOProject")
    # --- Excel Reads ---
    df = pd.read_excel(xlsx_path, sheet_name="Groups")
    df.columns = df.columns.str.strip()
    df['Group'] = df['Group'].astype(str).str.strip()
    df = df.set_index('Group')

    dfcp = pd.read_excel(xlsx_path, sheet_name="CpGC")
    dfcp.columns = dfcp.columns.str.strip()
    dfcp['Group'] = dfcp['Group'].astype(str).str.strip()
    dfcp = dfcp.set_index('Group')

    dfSVSH = pd.read_excel(xlsx_path, sheet_name="SVSH", index_col=0)
    dfSVSH.columns = dfSVSH.columns.str.strip()
    dfSVSH.index = dfSVSH.index.str.strip()

    # --- Max individual groups count (used as bounds and big-M for cuts) ---
    ni_max = 8
    #max total groups
    most=15
    least=1

    # --- Sets ---
    model.g = Set(initialize=dfcp.index.unique().tolist())
    model.i = Set(initialize=df.index.unique().tolist())

    A_props = ['Tci', 'Pci', 'Tbi', 'A0i', 'B0i', 'C0i', 'D0i']
    model.A = Set(initialize=A_props)

    X_props = ['Tb1i','Tm1i','δD1i','δP1i','δH1i','Vm1i']
    model.X = Set(initialize=X_props)

    model.B = Set(initialize=['acyclic','monocyclic'])

    # --- Parameters: Hukkerikar contributions (i,X) ---
    i_data = {(i, X): float(df.loc[i, X]) for i in df.index for X in X_props}
    model.ci = Param(model.i, model.X, initialize=i_data, default=0)

    # --- Parameters: CpGC contributions (g,A) ---
    g_data = {(g, A): float(dfcp.loc[g, A]) for g in dfcp.index for A in A_props}
    model.cg = Param(model.g, model.A, initialize=g_data, default=0)

    # --- Valency ---
    valency_dict = df['Valency 1'].to_dict()
    model.vi = Param(model.i, initialize=valency_dict, default=0)

    # --- Group types ---
    TYPE_COL = "Type of molecule (like aromatic and so on) ?"
    type_str_dict = df[TYPE_COL].to_dict()

    model.Ga = Set(initialize=[i for i, t in type_str_dict.items() if t == 'aromatic'])
    print("Ga groups:", ", ".join(sorted([str(i) for i in model.Ga])))
    model.Gc = Set(initialize=[i for i, t in type_str_dict.items() if t == 'cyclic'])
    model.Ggen = Set(initialize=[i for i, t in type_str_dict.items() if t == 'general'])
    model.Ggunsat = Set(initialize=[i for i, t in type_str_dict.items() if t == 'unsat'])

    # --- i -> g mapping (SVSH) ---
    ig_dict = {}
    for i in dfSVSH.index:
        for g in dfSVSH.columns:
            try:
                ig_dict[(i, g)] = int(dfSVSH.loc[i, g]) if pd.notna(dfSVSH.loc[i, g]) else 0
            except:
                ig_dict[(i, g)] = 0
    model.ig = Param(model.i, model.g, initialize=ig_dict, default=0)

    # --- Constants / Params ---
    model.Tb0 = Param(initialize=244.7889)   # K
    model.Tm0 = Param(initialize=144.0977)   # K
    model.Vm0 = Param(initialize=0.0123)     # m^3/kmol
    model.R0  = Param(initialize=3.3)        # MPa^0.5
    model.sigmacD = Param(initialize=15.7)
    model.sigmacP = Param(initialize=5.2)
    model.sigmacH = Param(initialize=5.8)
    model.T_avg   = Param(initialize=353)

    # Targets for standardized objective
    model.t_rho = Param(initialize=15.828)
    model.t_Cp  = Param(initialize=166.6)
    model.t_RED = Param(initialize=3.94)
    model.d_rho = Param(initialize=5)
    model.d_Cp  = Param(initialize=30)
    model.d_RED = Param(initialize=3)

    # Scaling params (not used in standardized form, kept for reference)
    model.rhomin = Param(initialize=0.1)
    model.rhomax = Param(initialize=36667)
    model.Cpmin  = Param(initialize=0.01)
    model.Cpmax  = Param(initialize=100)
    model.REDmin = Param(initialize=1e-7)
    model.REDmax = Param(initialize=10)
    model.eps    = Param(initialize=1e-12)

    # Objective weights
    model.w_RED = Param(initialize=1.0, mutable=True)
    model.w_rho = Param(initialize=1.0,  mutable=True)
    model.w_Cp  = Param(initialize=1.0,  mutable=True)

    # --- Variables (INTEGERS for counts) ---
    model.Tm       = Var(within=NonNegativeReals, bounds=(0.01, 313), initialize=5)
    model.Tb       = Var(within=NonNegativeReals, bounds=(393, 3000), initialize=500)
    model.rho      = Var(within=NonNegativeReals, bounds=(0.0001, 50000))
    model.RED      = Var(within=NonNegativeReals, bounds=(0.00001, None))
    model.sol1     = Var(within=Reals, bounds=(-16, 16))
    model.sol2     = Var(within=Reals, bounds=(-16, 16))
    model.sol3     = Var(within=Reals, bounds=(-16, 16))
    model.molarvol = Var(within=NonNegativeReals, bounds=(1e-6, 100000))

    # REMOVED binary expansion; use pure integers
    model.ni = Var(model.i, within=NonNegativeIntegers, bounds=(0, ni_max))
    model.ng = Var(model.g, within=NonNegativeIntegers, bounds=(0, 100))

    # Mode and selection binaries
    model.yb = Var(model.B, within=Binary)
    model.m  = Var(within=Reals, bounds=(-1, 2))
    model.ya = Var(within=Binary)
    model.yc = Var(within=Binary)

    # Cp intermediates
    model.T_b = Var(within=NonNegativeReals, initialize=350)
    model.T_c = Var(within=NonNegativeReals, initialize=500)
    model.P_c = Var(within=PositiveReals, initialize=1.0)
    model.T_avgr = Var(within=PositiveReals, initialize=0.7)
    model.W = Var(within=Reals, initialize=0.3)
    model.Cp0 = Var(within=NonNegativeReals, initialize=1)

    # Properties results
    model.Cp = Var(within=NonNegativeReals, bounds=(1e-6, 1000))
    model.fX = Var(model.X, within=Reals)

    # Big-M
    M_groups_val = ni_max * max(1, len(model.i))
    model.M_groups = Param(initialize=M_groups_val)

    # Indicators
    is_arom_dict = {i: 1 if i in list(model.Ga) else 0 for i in model.i}
    is_vgt2_arom_dict = {i: 1 if (i in list(model.Ga) and int(valency_dict.get(i, 0)) > 2) else 0 for i in model.i}
    model.is_arom      = Param(model.i, initialize=is_arom_dict, default=0)
    model.is_vgt2_arom = Param(model.i, initialize=is_vgt2_arom_dict, default=0)

    # Count expressions
    model.count_arom_vgt2    = Expression(expr=sum(model.ni[i] * model.is_vgt2_arom[i] for i in model.i))
    model.count_non_aromatic = Expression(expr=sum(model.ni[i] * (1 - model.is_arom[i]) for i in model.i))

    # Attachment flag
    model.z_attach_ok = Var(within=Binary)

    # --- Constraints ---
    def XGC_rule(model, X):
        return model.fX[X] == sum(model.ni[i] * model.ci[i, X] for i in model.i)
    model.XGC = Constraint(model.X, rule=XGC_rule)

    def numgini_rule(model, g):
        return model.ng[g] == sum(model.ni[i] * model.ig[i, g] for i in model.i)
    model.numgini = Constraint(model.g, rule=numgini_rule)

    # Link properties (recommend exp form to avoid ln domain issues)
    def link_properties_rule(model, X):
        if X == 'Tm1i':
            # Tm = Tm0 * exp(fX['Tm1i'])
            return model.fX[X] == exp(model.Tm/model.Tm0)
        elif X == 'Tb1i':
            # Tb = Tb0 * exp(fX['Tb1i'])
            return model.fX[X] == exp(model.Tb/model.Tb0)
        elif X == 'Vm1i':
            return model.molarvol == model.Vm0 + model.fX[X]
        elif X == 'δD1i':
            return model.sol1 == model.fX[X]
        elif X == 'δP1i':
            return model.sol2 == model.fX[X]
        elif X == 'δH1i':
            return model.sol3 == model.fX[X]
        else:
            return Constraint.Skip
    model.link_properties = Constraint(model.X, rule=link_properties_rule)

    # Density and RED
    model.rho_def = Constraint(expr=model.rho * model.molarvol == 1)
    model.RED_def = Constraint(expr=model.RED * model.R0 == sqrt(
        4 * (model.sol1 - model.sigmacD)**2 +
        (model.sol2 - model.sigmacP)**2 +
        (model.sol3 - model.sigmacH)**2
    ))

    # Cp correlation
    model.Tb_calc = Constraint(expr=model.T_b == 198.2 + sum(model.ng[g] * model.cg[g, 'Tbi'] for g in model.g))
    model.Tc_calc = Constraint(expr=model.T_c * (0.584 + 0.965 * sum(model.ng[g] * model.cg[g, 'Tci'] for g in model.g)) == model.T_b)
    model.Pc_calc = Constraint(expr=model.P_c == (0.113 + 0.0032 * sum(model.ng[g] for g in model.g) -
                                                sum(model.ng[g] * model.cg[g, 'Pci'] for g in model.g)))
    model.Tavgr_calc = Constraint(expr=model.T_avgr * model.T_c == model.T_avg)

    def W_calc_rule(model):
        alpha = (-5.97214 - log(model.P_c / 1.013 + (6.09648 * model.T_c) / model.T_b) +
                1.28862 * log(model.T_b / model.T_c) - 0.167347 * (model.T_b / model.T_c)**6)
        beta = (15.2518 - (15.6875 * model.T_c) / model.T_b - 13.4721 * log(model.T_b / model.T_c) +
                0.43577 * (model.T_b / model.T_c)**6)
        return model.W * beta == alpha
    model.W_calc = Constraint(rule=W_calc_rule)

    model.Cp0_calc = Constraint(expr=model.Cp0 == (
        sum(model.ng[g] * model.cg[g, 'A0i'] for g in model.g) - 37.0/93.0 +
        (sum(model.ng[g] * model.cg[g, 'B0i'] for g in model.g) + 0.21) * model.T_avg +
        (sum(model.ng[g] * model.cg[g, 'C0i'] for g in model.g) - 3.91e-4) * model.T_avg**2 +
        (sum(model.ng[g] * model.cg[g, 'D0i'] for g in model.g) + 2.06e-7) * model.T_avg**3
    ))

    model.Cp_constraint = Constraint(expr=model.Cp == (1/4.1868) * (
        model.Cp0 + 8.314 * (1.45 + (0.45/(1-model.T_avgr)) +
        0.25*model.W * (17.11 + 25.2*(((1 - model.T_avgr)**(1/3))/model.T_avgr) + (1.742 / (1-model.T_avgr)))))
    )

    # Mode selection and valency/bonding rules
    model.one_type = Constraint(expr=sum(model.yb[b] for b in model.B) == 1)
    model.m_rule   = Constraint(expr=model.m == 1 - model.yb['monocyclic'])

    def vi_rule(model):
        return 2*model.m == sum((2 - model.vi[i]) * model.ni[i] for i in model.i)
    model.vi_rule = Constraint(rule=vi_rule)

    model.monocyclic_mode_selection = Constraint(expr=model.ya + model.yc == model.yb['monocyclic'])
    model.aromatic_ring = Constraint(expr=sum(model.ni[i] for i in model.Ga) == 6 * model.ya)
    model.cyclic_only_if_cyclic_mode = Constraint(expr=sum(model.ni[i] for i in model.Gc) <= model.M_groups * model.yc)

    model.attach_ok_lower = Constraint(expr=model.count_arom_vgt2 >= model.z_attach_ok)
    model.attach_ok_upper = Constraint(expr=model.count_arom_vgt2 <= model.M_groups * model.z_attach_ok)
    model.non_aromatic_allowed_in_aromatic_mode = Constraint(
        expr=model.count_non_aromatic <= model.M_groups * (1 - model.ya) + model.M_groups * model.z_attach_ok
    )

    # At least n groups
    model.at_least = Constraint(expr=sum(model.ni[i] for i in model.i) >= least)

    # At least most n groups
    model.at_most = Constraint(expr=sum(model.ni[i] for i in model.i) <= most)
    # Min 5 cyclic groups if cyclic non-aromatic
    model.min_cyclic_groups = Constraint(expr=sum(model.ni[i] for i in model.Gc) >= 5 * model.yc)

    # ---- Exclusion (nogood) cut: robust build over intersection with model.i ----
    '''
    excl_combo = {
        'CH3': 1,
        'NH2 except as above': 1,
        '-O-': 2,
    }

    # Intersect requested names with model.i (robust to mismatches)

    excl_indices = [i for i in excl_combo.keys() if i in list(model.i)]
    model.EXCL = Set(initialize=excl_indices)

    # If you want to fail fast when a name is missing, uncomment:
    # missing = [i for i in excl_combo.keys() if i not in list(model.i)]
    # if missing:
    #     raise ValueError(f"Nogood cut: groups not found in model.i: {missing}. Available: {list(model.i)}")

    # Target counts as a Param over EXCL
    model.target_excl = Param(model.EXCL, initialize={i: excl_combo[i] for i in excl_indices}, default=0)

    # Binary indicators: diff[i] = 1 if ni[i] differs from target_excl[i], otherwise 0
    model.diff = Var(model.EXCL, within=Binary)

    # Safe Big-M: use ni_max
    M_cut = ni_max

    def diff_upper_rule(model, i):
        # ni[i] - target[i] <= M * diff[i]
        return model.ni[i] - model.target_excl[i] <= M_cut * model.diff[i]
    model.diff_upper = Constraint(model.EXCL, rule=diff_upper_rule)

    def diff_lower_rule(model, i):
        # target[i] - ni[i] <= M * diff[i]
        return model.target_excl[i] - model.ni[i] <= M_cut * model.diff[i]
    model.diff_lower = Constraint(model.EXCL, rule=diff_lower_rule)

    def exclude_exact_match_rule(model):
        # If EXCL is empty, skip the cut to avoid a 0 >= 1 boolean
        if len(model.EXCL) == 0:
            return Constraint.Skip
        return sum(model.diff[i] for i in model.EXCL) >= 1
    model.exclude_exact_match = Constraint(rule=exclude_exact_match_rule)
    '''
    # --- Standardized objective (targets and spreads) ---
    model.rho_norm = Expression(expr=(model.rho - model.t_rho) / (model.d_rho))
    model.Cp_norm  = Expression(expr=(model.Cp  - model.t_Cp)  / (model.d_Cp))
    model.RED_norm = Expression(expr=(model.RED - model.t_RED) / (model.d_RED))

    def objective_rule(model):
        return model.w_RED * model.RED_norm + model.w_rho * (1 - model.rho_norm) + model.w_Cp * model.Cp_norm
    model.obj = Objective(rule=objective_rule, sense=minimize)

    return model


# ------------------------------------------------------------
# 2) Runner: set weights, solve, extract results
# ------------------------------------------------------------
def run_with_weights(base_model: ConcreteModel, weights, solver_name="gams", gams_solver="DICOPT"):
    # Clone model to avoid cross-run interference
    m = base_model.clone()

    # Set weights
    w_red, w_rho, w_cp = weights
    m.w_RED.set_value(float(w_red))
    m.w_rho.set_value(float(w_rho))
    m.w_Cp.set_value(float(w_cp))

    # Solve
    solver = SolverFactory(solver_name)
    if solver_name == "gams":
        solver.options['solver'] = gams_solver
    results = solver.solve(m, tee=False)

    # Feasibility check
    term = results.solver.termination_condition
    feasible = term in {TerminationCondition.optimal, TerminationCondition.feasible}
    status = str(results.solver.status)
    termination = str(term)

    # Helper
    def v(x):
        try:
            return float(value(x))
        except:
            return None

    # Mode string
    chosen_type = None
    for b in m.B:
        try:
            if v(m.yb[b]) is not None and v(m.yb[b]) > 0.5:
                chosen_type = b
                break
        except:
            pass

    # Compose outputs
    metrics = {
        "w_RED": w_red, "w_rho": w_rho, "w_Cp": w_cp,
        "obj": v(m.obj),
        "RED_norm": v(m.RED_norm),
        "rho_norm": v(m.rho_norm),
        "Cp_norm":  v(m.Cp_norm),
        "Cp":  v(m.Cp),
        "rho": v(m.rho),
        "RED": v(m.RED),
        "Tm":  v(m.Tm),
        "Tb":  v(m.Tb),
        "Vm":  v(m.molarvol),
        "T_b": v(m.T_b),
        "T_c": v(m.T_c),
        "P_c": v(m.P_c),
        "T_avgr": v(m.T_avgr),
        "W": v(m.W),
        "Cp0": v(m.Cp0),
        "sol1": v(m.sol1),
        "sol2": v(m.sol2),
        "sol3": v(m.sol3),
        "mode": chosen_type,
        "m_param": v(m.m),
        "solver_status": status,
        "termination": termination,
        "feasible": feasible,
    }

    # Compositions (integers)
    ni = {str(i): int(round(v(m.ni[i]))) if v(m.ni[i]) is not None else None for i in m.i}
    ng = {str(g): int(round(v(m.ng[g]))) if v(m.ng[g]) is not None else None for g in m.g}

    return metrics, ni, ng


# ------------------------------------------------------------
# 3) Loop over weight vectors and export to Excel
# ------------------------------------------------------------
def run_batch_and_export(xlsx_path: str, vectors, out_xlsx="results_weights.xlsx",
                         solver_name="gams", gams_solver="DICOPT"):
    base_model = build_model(xlsx_path)

    metrics_rows = []
    ni_tables = []
    ng_tables = []

    for idx, w in enumerate(vectors, start=1):
        metrics, ni, ng = run_with_weights(base_model, w, solver_name, gams_solver)
        metrics["run_id"] = idx
        metrics_rows.append(metrics)

        # Build ni/ng DataFrames for this run
        ni_df = pd.DataFrame({"group": list(ni.keys()), f"count_run_{idx}": list(ni.values())})
        ng_df = pd.DataFrame({"group": list(ng.keys()), f"count_run_{idx}": list(ng.values())})
        ni_tables.append(ni_df)
        ng_tables.append(ng_df)

    # Metrics summary
    metrics_df = pd.DataFrame(metrics_rows)

    # Merge ni/ng tables across runs (outer-join on group)
    ni_merged = None
    for df in ni_tables:
        ni_merged = df if ni_merged is None else ni_merged.merge(df, on="group", how="outer")
    ng_merged = None
    for df in ng_tables:
        ng_merged = df if ng_merged is None else ng_merged.merge(df, on="group", how="outer")

    # Export to Excel with multiple sheets
    with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
        metrics_df.to_excel(writer, sheet_name="metrics", index=False)
        if ni_merged is not None:
            ni_merged.to_excel(writer, sheet_name="ni", index=False)
        if ng_merged is not None:
            ng_merged.to_excel(writer, sheet_name="ng", index=False)

    return metrics_df, ni_merged, ng_merged


# ------------------------------------------------------------
# 4) Example usage
# ------------------------------------------------------------
if __name__ == "__main__":
    xlsx_path = r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx"

    vectors = [
        [1.0, 1.0, 1.0],
        [1.5, 1.5, 0.5],
        [0.5, 1.5, 1.5],
        [1.5, 0.5, 1.5],
        [0.5, 0.5, 1.5],
        [1.5, 0.5, 0.5],
    ]

    metrics_df, ni_df, ng_df = run_batch_and_export(
        xlsx_path=xlsx_path,
        vectors=vectors,
        out_xlsx="results_weights.xlsx",
        solver_name="gams",
        gams_solver="DICOPT"  # try 'antigone' if DICOPT is infeasible
    )

    print("Saved results to results_weights.xlsx")
    print(metrics_df.head())
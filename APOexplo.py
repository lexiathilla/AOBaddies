#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Standalone feasibility sampler to compute mean/std for optimization:
- Reads Groups, CpGC, SVSH sheets from your Excel file
- Builds a feasibility-only Pyomo model
- Randomly samples ni counts to generate feasible molecules
- Computes min–max normalized metrics and z-scores (mean/std from sample)
- Saves stats to JSON/CSV for downstream optimization

Usage:
  python sample_stats.py --xlsx "C:\\path\\to\\GCs.xlsx" --samples 100 --seed 123 --ni_max 3
"""

from pyomo.environ import (
    ConcreteModel, Set, Param, Var, RangeSet,
    NonNegativeReals, PositiveReals, Reals,
    Expression, Constraint, Objective, minimize,
    log, exp, sqrt, value, SolverFactory
)
from pyomo.opt import TerminationCondition
import pandas as pd
import math
import random
import argparse
import json
import sys


# ---------- Data reading ----------
def read_sheets(xlsx_path: str):
    # Groups
    df = pd.read_excel(xlsx_path, sheet_name="Groups")
    df.columns = df.columns.str.strip()
    if "Group" not in df.columns:
        raise ValueError("Groups sheet must contain a 'Group' column.")
    df["Group"] = df["Group"].astype(str).str.strip()
    df = df.set_index("Group")

    # CpGC
    dfcp = pd.read_excel(xlsx_path, sheet_name="CpGC")
    dfcp.columns = dfcp.columns.str.strip()
    if "Group" not in dfcp.columns:
        raise ValueError("CpGC sheet must contain a 'Group' column.")
    dfcp["Group"] = dfcp["Group"].astype(str).str.strip()
    dfcp = dfcp.set_index("Group")

    # SVSH (i->g mapping)
    dfSVSH = pd.read_excel(xlsx_path, sheet_name="SVSH", index_col=0)
    dfSVSH.columns = dfSVSH.columns.str.strip()
    dfSVSH.index = dfSVSH.index.astype(str).str.strip()

    # Clean placeholders
    for d in (df, dfcp, dfSVSH):
        d.replace("****", pd.NA, inplace=True)

    return df, dfcp, dfSVSH


# ---------- Build feasibility model ----------
def build_model(df: pd.DataFrame, dfcp: pd.DataFrame, dfSVSH: pd.DataFrame,
                ni_max: int = 3,
                tb0: float = 244.7889, tm0: float = 144.0977,
                vm0: float = 0.0123, r0: float = 4.7,
                sigD: float = 15.6, sigP: float = 5.2, sigH: float = 5.8,
                t_avg: float = 353.0,
                rhomin: float = 0.1, rhomax: float = 36667.0,
                cpmin: float = 0.01, cpmax: float = 100.0,
                redmin: float = 1e-7, redmax: float = 10.0):
    model = ConcreteModel("FeasibilitySampler")

    # Sets
    model.i = Set(initialize=df.index.unique().tolist())
    model.g = Set(initialize=dfcp.index.unique().tolist())

    X_props = ["Tb1i", "Tm1i", "δD1i", "δP1i", "δH1i", "Vm1i"]
    A_props = ["Tci", "Pci", "Tbi", "A0i", "B0i", "C0i", "D0i"]
    for col in X_props:
        if col not in df.columns:
            raise ValueError(f"Groups sheet missing column: {col}")
    for col in A_props:
        if col not in dfcp.columns:
            raise ValueError(f"CpGC sheet missing column: {col}")
    model.X = Set(initialize=X_props)
    model.A = Set(initialize=A_props)

    # Params: Hukkerikar (i,X)
    df[X_props] = df[X_props].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    i_data = {(i, X): float(df.loc[i, X]) for i in df.index for X in X_props}
    model.ci = Param(model.i, model.X, initialize=i_data, default=0.0)

    # Params: CpGC (g,A)
    dfcp[A_props] = dfcp[A_props].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    g_data = {(g, A): float(dfcp.loc[g, A]) for g in dfcp.index for A in A_props}
    model.cg = Param(model.g, model.A, initialize=g_data, default=0.0)

    # SVSH i->g map (use only CpGC-group columns)
    svsh_cols = [c for c in dfSVSH.columns if c in model.g]
    df_svsh_num = dfSVSH[svsh_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    ig_dict = {}
    for i_name in dfSVSH.index:
        if i_name not in model.i:
            continue
        for g_name in svsh_cols:
            ig_dict[(i_name, g_name)] = int(df_svsh_num.loc[i_name, g_name])
    model.ig = Param(model.i, model.g, initialize=ig_dict, default=0)

    # Constants
    model.Tb0 = Param(initialize=tb0)   # K
    model.Tm0 = Param(initialize=tm0)   # K
    model.Vm0 = Param(initialize=vm0)   # m^3/kmol
    model.R0  = Param(initialize=r0)    # MPa^0.5
    model.sigmacD = Param(initialize=sigD)
    model.sigmacP = Param(initialize=sigP)
    model.sigmacH = Param(initialize=sigH)
    model.T_avg   = Param(initialize=t_avg)

    # Normalization params (used to compute normalized metrics)
    model.rhomin = Param(initialize=rhomin)
    model.rhomax = Param(initialize=rhomax)
    model.Cpmin  = Param(initialize=cpmin)
    model.Cpmax  = Param(initialize=cpmax)
    model.REDmin = Param(initialize=redmin)
    model.REDmax = Param(initialize=redmax)
    model.eps    = Param(initialize=1e-12)

    # Variables
    model.ni = Var(model.i, within=NonNegativeReals, bounds=(0, ni_max))
    model.ng = Var(model.g, within=NonNegativeReals, bounds=(0, None))
    model.fX = Var(model.X, within=Reals)

    model.Tm       = Var(within=NonNegativeReals, bounds=(1e-2, 313))
    model.Tb       = Var(within=NonNegativeReals, bounds=(1e-2, 3000))   # relaxed lower bound to improve feasibility
    model.molarvol = Var(within=NonNegativeReals, bounds=(1e-6, 1e5))
    model.rho      = Var(within=NonNegativeReals, bounds=(1e-6, None))
    model.RED      = Var(within=NonNegativeReals, bounds=(1e-6, None))
    model.sol1     = Var(within=Reals, bounds=(-16, 16))
    model.sol2     = Var(within=Reals, bounds=(-16, 16))
    model.sol3     = Var(within=Reals, bounds=(-16, 16))

    # Cp intermediates
    model.T_b   = Var(within=NonNegativeReals, initialize=350)          # K
    model.T_c   = Var(within=NonNegativeReals, initialize=500)          # K
    model.P_c   = Var(within=PositiveReals, initialize=1.0)             # MPa
    model.T_avgr= Var(within=PositiveReals, bounds=(1e-3, 0.95), initialize=0.7)
    model.W     = Var(within=Reals, initialize=0.3)
    model.Cp0   = Var(within=NonNegativeReals, initialize=1)
    model.Cp    = Var(within=NonNegativeReals, bounds=(1e-6, 1000))

    # Constraints
    def XGC_rule(m, X):
        return m.fX[X] == sum(m.ni[i] * m.ci[i, X] for i in m.i)
    model.XGC = Constraint(model.X, rule=XGC_rule)

    def numgini_rule(m, g):
        return m.ng[g] == sum(m.ni[i] * m.ig[i, g] for i in m.i)
    model.numgini = Constraint(model.g, rule=numgini_rule)

    # Hukkerikar property links (exp form to avoid ln-domain issues)
    def link_properties_rule(m, X):
        if X == "Tm1i":
            return m.Tm == m.Tm0 * exp(m.fX[X])
        elif X == "Tb1i":
            return m.Tb == m.Tb0 * exp(m.fX[X])
        elif X == "Vm1i":
            return m.molarvol == m.Vm0 + m.fX[X]
        elif X == "δD1i":
            return m.sol1 == m.fX[X]
        elif X == "δP1i":
            return m.sol2 == m.fX[X]
        elif X == "δH1i":
            return m.sol3 == m.fX[X]
        return Constraint.Skip
    model.link_properties = Constraint(model.X, rule=link_properties_rule)

    # Density and RED
    model.rho_def = Constraint(expr=model.rho * model.molarvol == 1.0)
    model.RED_def = Constraint(expr=model.RED * model.R0 == sqrt(
        4.0 * (model.sol1 - model.sigmacD)**2 +
        (model.sol2 - model.sigmacP)**2 +
        (model.sol3 - model.sigmacH)**2
    ))

    # Cp block (Sahinidis/Rowlinson)
    model.Tb_calc = Constraint(expr=model.T_b == 198.2 + sum(model.ng[g] * model.cg[g, "Tbi"] for g in model.g))
    model.Tc_calc = Constraint(expr=model.T_c * (0.584 + 0.965 * sum(model.ng[g] * model.cg[g, "Tci"] for g in model.g)) == model.T_b)
    model.Pc_calc = Constraint(expr=model.P_c == (0.113 + 0.0032 * sum(model.ng[g] for g in model.g) -
                                                 sum(model.ng[g] * model.cg[g, "Pci"] for g in model.g)))
    model.Tavgr_calc = Constraint(expr=model.T_avgr * model.T_c == model.T_avg)

    def W_calc_rule(m):
        alpha = (-5.97214 - log(m.P_c / 1.013 + (6.09648 * m.T_c) / m.T_b) +
                 1.28862 * log(m.T_b / m.T_c) - 0.167347 * (m.T_b / m.T_c)**6)
        beta = (15.2518 - (15.6875 * m.T_c) / m.T_b - 13.4721 * log(m.T_b / m.T_c) +
                0.43577 * (m.T_b / m.T_c)**6)
        return m.W * beta == alpha
    model.W_calc = Constraint(rule=W_calc_rule)

    model.Cp0_calc = Constraint(expr=model.Cp0 == (
        sum(model.ng[g] * model.cg[g, "A0i"] for g in model.g) - 37.0/93.0 +
        (sum(model.ng[g] * model.cg[g, "B0i"] for g in model.g) + 0.21) * model.T_avg +
        (sum(model.ng[g] * model.cg[g, "C0i"] for g in model.g) - 3.91e-4) * model.T_avg**2 +
        (sum(model.ng[g] * model.cg[g, "D0i"] for g in model.g) + 2.06e-7) * model.T_avg**3
    ))

    model.Cp_constraint = Constraint(expr=model.Cp == (1/4.1868) * (
        model.Cp0 + 8.314 * (1.45 + (0.45 / (1 - model.T_avgr)) +
        0.25 * model.W * (17.11 + 25.2 * (((1 - model.T_avgr)**(1/3)) / model.T_avgr) + (1.742 / (1 - model.T_avgr)))))
    ))

    # At least one group
    model.at_least_one = Constraint(expr=sum(model.ni[i] for i in model.i) >= 1)

    # Normalized Expressions (min–max)
    model.rho_norm = Expression(expr=(model.rho - model.rhomin) / (model.rhomax - model.rhomin + model.eps))
    model.Cp_norm  = Expression(expr=(model.Cp  - model.Cpmin)  / (model.Cpmax  - model.Cpmin  + model.eps))
    model.RED_norm = Expression(expr=(log(model.RED) - log(model.REDmin)) / (log(model.REDmax) - log(model.REDmin) + model.eps))

    # Feasibility objective (constant)
    model.feas_obj = Objective(expr=0.0)

    return model


# ---------- Randomize counts ----------
def randomize_counts(model, ni_max: int, rng: random.Random):
    total = 0
    for i in model.i:
        val = rng.randint(0, ni_max)
        model.ni[i].fix(val)
        total += val
    if total == 0:
        # ensure at least one group
        i_pick = rng.choice(list(model.i.value))
        model.ni[i_pick].fix(1)


# ---------- Solve one sample ----------
def solve_one_sample(base_model, ni_max, rng, solver_name="gams", gams_solver="antigone"):
    m = base_model.clone()
    randomize_counts(m, ni_max, rng)

    solver = SolverFactory(solver_name)
    if not solver.available():
        if SolverFactory("ipopt").available():
            solver_name = "ipopt"
            solver = SolverFactory("ipopt")
        else:
            raise RuntimeError("No available solver found (gams or ipopt).")
    if solver_name == "gams":
        solver.options["solver"] = gams_solver

    results = solver.solve(m, tee=False)
    term = results.solver.termination_condition
    feasible = term in {TerminationCondition.optimal, TerminationCondition.feasible}
    if not feasible:
        return None

    try:
        rec = {
            "rho": float(value(m.rho)),
            "Cp": float(value(m.Cp)),
            "RED": float(value(m.RED)),
            "rho_norm": float(value(m.rho_norm)),
            "Cp_norm": float(value(m.Cp_norm)),
            "RED_norm": float(value(m.RED_norm)),
            "Tm": float(value(m.Tm)),
            "Tb": float(value(m.Tb)),
            "Vm": float(value(m.molarvol)),
            "ni": {str(i): float(value(m.ni[i])) for i in m.i},
            "ng": {str(g): float(value(m.ng[g])) for g in m.g},
        }
    except Exception:
        return None

    return rec


# ---------- Sample N feasible molecules and compute stats ----------
def sample_and_stats(base_model, n_samples=100, ni_max=3, seed=123,
                     solver_name="gams", gams_solver="antigone", tries_per_sample=8):
    rng = random.Random(seed)
    records = []
    for s in range(n_samples):
        sample = None
        for _ in range(tries_per_sample):
            sample = solve_one_sample(base_model, ni_max, rng, solver_name, gams_solver)
            if sample is not None:
                break
        if sample is None:
            continue
        records.append(sample)

    if not records:
        raise RuntimeError("No feasible samples generated. Relax bounds/constraints or increase tries_per_sample.")

    df = pd.DataFrame([{
        "rho": r["rho"], "Cp": r["Cp"], "RED": r["RED"],
        "rho_norm": r["rho_norm"], "Cp_norm": r["Cp_norm"], "RED_norm": r["RED_norm"],
        "Tm": r["Tm"], "Tb": r["Tb"], "Vm": r["Vm"],
        "ni": r["ni"], "ng": r["ng"]
    } for r in records])

    # Means and stds for normalized metrics (for use in optimization)
    stats = {}
    for col in ["rho_norm", "Cp_norm", "RED_norm"]:
        mu = float(df[col].mean())
        sd = float(df[col].std(ddof=1)) if df.shape[0] > 1 else 0.0
        stats[col] = {"mean": mu, "std": sd}

    # Also provide raw metric stats (optional, in case you prefer standardizing raw)
    for col in ["rho", "Cp", "RED"]:
        mu = float(df[col].mean())
        sd = float(df[col].std(ddof=1)) if df.shape[0] > 1 else 0.0
        stats[col] = {"mean": mu, "std": sd}

    return df, stats


# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser(description="Generate feasible samples and compute mean/std for optimization")
    parser.add_argument("--xlsx", type=str, required=True, help="Path to GCs.xlsx")
    parser.add_argument("--samples", type=int, default=100, help="Number of feasible molecules to sample")
    parser.add_argument("--seed", type=int, default=123, help="Random seed")
    parser.add_argument("--ni_max", type=int, default=3, help="Max count per Hukkerikar group")
    parser.add_argument("--solver", type=str, default="gams", choices=["gams", "ipopt"], help="Solver to use")
    parser.add_argument("--gams_solver", type=str, default="antigone", help="GAMS NLP/MINLP solver (e.g., antigone, conopt)")
    # Normalization params used during sampling
    parser.add_argument("--rhomin", type=float, default=0.1)
    parser.add_argument("--rhomax", type=float, default=36667.0)
    parser.add_argument("--cpmin", type=float, default=0.01)
    parser.add_argument("--cpmax", type=float, default=100.0)
    parser.add_argument("--redmin", type=float, default=1e-7)
    parser.add_argument("--redmax", type=float, default=10.0)
    parser.add_argument("--out_csv", type=str, default="feasible_samples.csv", help="Output CSV with samples")
    parser.add_argument("--out_json", type=str, default="norm_stats.json", help="Output JSON with mean/std")
    args = parser.parse_args()

    try:
        df, dfcp, dfSVSH = read_sheets(args.xlsx)
    except Exception as e:
        print(f"Error reading sheets: {e}", file=sys.stderr)
        sys.exit(1)

    model = build_model(df, dfcp, dfSVSH,
                        ni_max=args.ni_max,
                        rhomin=args.rhomin, rhomax=args.rhomax,
                        cpmin=args.cpmin,   cpmax=args.cpmax,
                        redmin=args.redmin, redmax=args.redmax)

    df_samples, stats = sample_and_stats(model,
                                         n_samples=args.samples,
                                         ni_max=args.ni_max,
                                         seed=args.seed,
                                         solver_name=args.solver,
                                         gams_solver=args.gams_solver)

    # Save artifacts
    df_samples.to_csv(args.out_csv, index=False)
    with open(args.out_json, "w") as f:
        json.dump(stats, f, indent=2)

    print(f"\nSaved {len(df_samples)} feasible samples to: {args.out_csv}")
    print(f"Saved normalization stats (mean/std) to: {args.out_json}")
    print("\nNormalization stats (for optimization):")
    for k, v in stats.items():
        print(f"  {k}: mean = {v['mean']:.6f}, std = {v['std']:.6f}")


if __name__ == "__main__":
    main()
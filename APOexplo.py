from pyomo.environ import *
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random

inputGCs = r"C:\Users\Alexia\OneDrive - Imperial College London\AAYEAR4\APO1\GCs.xlsx"

# ============================================================================
# EXCEL READS
# ============================================================================
df = pd.read_excel(inputGCs, sheet_name="Groups")
df.columns = df.columns.str.strip()
df['Group'] = df['Group'].astype(str).str.strip()
df = df.set_index('Group')

dfcp = pd.read_excel(inputGCs, sheet_name="CpGC")
dfcp.columns = dfcp.columns.str.strip()
dfcp['Group'] = dfcp['Group'].astype(str).str.strip()
dfcp = dfcp.set_index('Group')

dfSVSH = pd.read_excel(inputGCs, sheet_name="SVSH", index_col=0)
dfSVSH.columns = dfSVSH.columns.str.strip()
dfSVSH.index = dfSVSH.index.str.strip()

# ============================================================================
# PARAMETERS
# ============================================================================
ni_max = 8
X_props = ['Tb1i', 'Tm1i', 'δD1i', 'δP1i', 'δH1i', 'Vm1i']
A_props = ['Tci', 'Tbi', 'Pci', 'A0i', 'B0i', 'C0i', 'D0i']

# Physical constants
Tb0 = 244.7889  # K
Tm0 = 144.0977  # K
Vm0 = 20.7339   # m^3/kmol
R0 = 4.7        # MPa^0.5
sigmacD = 15.6  # MPa^0.5
sigmacP = 5.2   # MPa^0.5
sigmacH = 5.8   # MPa^0.5
T_avg = 353     # K

# Extract data
ci = {(i, X): float(df.loc[i, X]) for i in df.index for X in X_props}
cg = {(g, A): float(dfcp.loc[g, A]) for g in dfcp.index for A in A_props}
valency_dict = df['Valency 1'].to_dict()

# SVSH mapping
ig_dict = {}
for i in dfSVSH.index:
    for g in dfSVSH.columns:
        try:
            ig_dict[(i, g)] = int(dfSVSH.loc[i, g]) if pd.notna(dfSVSH.loc[i, g]) else 0
        except:
            ig_dict[(i, g)] = 0

# Group types
TYPE_COL = "Type of molecule (like aromatic and so on) ?"
type_str_dict = df[TYPE_COL].to_dict()
Ga = [i for i, t in type_str_dict.items() if t == 'aromatic']
Gc = [i for i, t in type_str_dict.items() if t == 'cyclic']
Ggen = [i for i, t in type_str_dict.items() if t == 'general']

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def calculate_properties(ni_dict, verbose=False):
    """Calculate all properties for a given group count configuration"""
    try:
        # Calculate fX values (Hukkerikar contributions)
        fX = {}
        for X in X_props:
            fX[X] = sum(ni_dict.get(i, 0) * ci[(i, X)] for i in df.index)
        
        # Calculate ng values (Sahinidis groups from Hukkerikar groups)
        ng = {}
        for g in dfcp.index:
            ng[g] = sum(ni_dict.get(i, 0) * ig_dict.get((i, g), 0) for i in df.index)
        
        # Link properties (Hukkerikar)
        Tm = np.exp(fX['Tm1i'] / Tm0)
        Tb = np.exp(fX['Tb1i'] / Tb0)
        molarvol = fX['Vm1i'] + Vm0
        sol1 = fX['δD1i']
        sol2 = fX['δP1i']
        sol3 = fX['δH1i']
        
        # Density
        rho = 1 / molarvol if molarvol > 0 else np.nan
        
        # RED
        RED = np.sqrt(4 * (sol1 - sigmacD)**2 + 
                      (sol2 - sigmacP)**2 + 
                      (sol3 - sigmacH)**2) / R0
        
        # Cp calculations (Sahinidis)
        T_b = 198.2 + sum(ng[g] * cg[(g, 'Tbi')] for g in dfcp.index)
        
        denom_tc = 0.584 + 0.965 * sum(ng[g] * cg[(g, 'Tci')] for g in dfcp.index)
        T_c = T_b / denom_tc if denom_tc > 0 else np.nan
        
        P_c = (0.113 + 0.0032 * sum(ng.values()) - 
               sum(ng[g] * cg[(g, 'Pci')] for g in dfcp.index))
        
        if T_c > 0 and T_c > T_avg:
            T_avgr = T_avg / T_c
        else:
            T_avgr = np.nan
        
        # Acentric factor
        if T_c > 0 and T_b > 0 and P_c > 0 and T_b != T_c:
            alpha = (-5.97214 - np.log(P_c / 1.013 + (6.09648 * T_c) / T_b) + 
                     1.28862 * np.log(T_b / T_c) - 0.167347 * (T_b / T_c)**6)
            beta = (15.2518 - (15.6875 * T_c) / T_b - 13.4721 * np.log(T_b / T_c) + 
                    0.43577 * (T_b / T_c)**6)
            W = alpha / beta if beta != 0 else np.nan
        else:
            W = np.nan
        
        # Cp0
        Cp0 = (sum(ng[g] * cg[(g, 'A0i')] for g in dfcp.index) - 37.0/93.0 +
               (sum(ng[g] * cg[(g, 'B0i')] for g in dfcp.index) + 0.21) * T_avg +
               (sum(ng[g] * cg[(g, 'C0i')] for g in dfcp.index) - 3.91e-4) * T_avg**2 +
               (sum(ng[g] * cg[(g, 'D0i')] for g in dfcp.index) + 2.06e-7) * T_avg**3)
        
        # Final Cp (Rowlinson)
        if not np.isnan(T_avgr) and T_avgr < 1 and not np.isnan(W):
            term1 = 1.45
            term2 = 0.45 / (1 - T_avgr)
            term3 = 0.25 * W * (17.11 + 25.2 * ((1 - T_avgr)**(1/3) / T_avgr) + 
                                1.742 / (1 - T_avgr))
            Cp = (1/4.1868) * (Cp0 + 8.314 * (term1 + term2 + term3))
        else:
            Cp = np.nan
        
        results = {
            'ni_dict': ni_dict.copy(),
            'ng_dict': ng.copy(),
            'Tm': Tm,
            'Tb': Tb,
            'rho': rho,
            'RED': RED,
            'Cp': Cp,
            'molarvol': molarvol,
            'sol1': sol1,
            'sol2': sol2,
            'sol3': sol3,
            'T_b': T_b,
            'T_c': T_c,
            'P_c': P_c,
            'T_avgr': T_avgr,
            'W': W,
            'Cp0': Cp0,
            'total_groups': sum(ni_dict.values()),
            'num_group_types': sum(1 for v in ni_dict.values() if v > 0)
        }
        
        return results
        
    except Exception as e:
        if verbose:
            print(f"Error in calculation: {e}")
        return None

def generate_random_molecule(ni_max=8, aromatic=False, cyclic=False):
    """Generate random group counts"""
    ni_dict = {}
    available_groups = list(df.index)
    
    if aromatic:
        # Must have exactly 6 aromatic groups
        for i in Ga:
            ni_dict[i] = 6 / len(Ga) if len(Ga) > 0 else 0
    elif cyclic:
        # Random cyclic groups
        available_groups = Gc
    else:
        # Acyclic - general groups
        available_groups = Ggen
    
    # Add random groups
    num_group_types = random.randint(2, min(6, len(available_groups)))
    selected_groups = random.sample(available_groups, num_group_types)
    
    for i in selected_groups:
        ni_dict[i] = random.randint(1, ni_max)
    
    return ni_dict

# ============================================================================
# SAMPLING
# ============================================================================

print("\n" + "="*80)
print("PROPERTY SAMPLING FOR SCALING ANALYSIS")
print("="*80)

num_samples = 100  # Number of random molecules to generate
results_list = []

print(f"\nGenerating {num_samples} random molecules...")
print("Progress: ", end="", flush=True)

for i in range(num_samples):
    if i % 10 == 0:
        print(f"{i}...", end="", flush=True)
    
    # Generate random molecule configuration
    molecule_type = random.choice(['acyclic', 'aromatic', 'cyclic'])
    
    if molecule_type == 'aromatic':
        ni_dict = generate_random_molecule(ni_max=ni_max, aromatic=True)
    elif molecule_type == 'cyclic':
        ni_dict = generate_random_molecule(ni_max=ni_max, cyclic=True)
    else:
        ni_dict = generate_random_molecule(ni_max=ni_max)
    
    # Calculate properties
    result = calculate_properties(ni_dict)
    
    if result is not None:
        result['molecule_type'] = molecule_type
        result['sample_id'] = i
        results_list.append(result)

print(f"\n\nSuccessfully calculated properties for {len(results_list)} molecules")

# ============================================================================
# ANALYZE RESULTS
# ============================================================================

# Convert to DataFrame for analysis
df_results = pd.DataFrame(results_list)

# Filter out NaN/Inf values
valid_results = df_results[
    np.isfinite(df_results['Cp']) & 
    np.isfinite(df_results['rho']) & 
    np.isfinite(df_results['RED']) &
    (df_results['Cp'] > 0) &
    (df_results['rho'] > 0) &
    (df_results['RED'] > 0)
]

print(f"Valid results (no NaN/Inf): {len(valid_results)}/{len(results_list)}")

if len(valid_results) > 0:
    print("\n" + "="*80)
    print("PROPERTY STATISTICS")
    print("="*80)
    
    properties = ['Cp', 'rho', 'RED', 'Tm', 'Tb', 'molarvol', 'sol1', 'sol2', 'sol3']
    
    for prop in properties:
        if prop in valid_results.columns:
            values = valid_results[prop]
            print(f"\n{prop}:")
            print(f"  Min:    {values.min():.6f}")
            print(f"  Max:    {values.max():.6f}")
            print(f"  Mean:   {values.mean():.6f}")
            print(f"  Median: {values.median():.6f}")
            print(f"  Std:    {values.std():.6f}")
    
    print("\n" + "="*80)
    print("RECOMMENDED SCALING BOUNDS")
    print("="*80)
    
    # Add safety margins
    print(f"\nHeat Capacity (Cp):")
    print(f"  Cpmin: {valid_results['Cp'].min() * 0.9:.6f}")
    print(f"  Cpmax: {valid_results['Cp'].max() * 1.1:.6f}")
    
    print(f"\nDensity (rho):")
    print(f"  rhomin: {valid_results['rho'].min() * 0.9:.2f}")
    print(f"  rhomax: {valid_results['rho'].max() * 1.1:.2f}")
    
    print(f"\nRED:")
    print(f"  REDmin: {valid_results['RED'].min() * 0.9:.6f}")
    print(f"  REDmax: {valid_results['RED'].max() * 1.1:.6f}")
    
    print(f"\nMolar Volume:")
    print(f"  Vm_min: {valid_results['molarvol'].min() * 0.9:.6f}")
    print(f"  Vm_max: {valid_results['molarvol'].max() * 1.1:.6f}")
    
    # ========================================================================
    # VISUALIZATION
    # ========================================================================
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle('Property Distributions from Random Sampling', fontsize=16)
    
    plot_props = ['Cp', 'rho', 'RED', 'Tm', 'Tb', 'molarvol', 'sol1', 'sol2', 'sol3']
    
    for idx, prop in enumerate(plot_props):
        ax = axes[idx // 3, idx % 3]
        if prop in valid_results.columns:
            valid_results[prop].hist(bins=30, ax=ax, edgecolor='black')
            ax.set_xlabel(prop)
            ax.set_ylabel('Frequency')
            ax.set_title(f'{prop} Distribution')
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('property_distributions.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved histogram to 'property_distributions.png'")
    
    # Scatter plots for objective components
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle('Key Property Relationships', fontsize=16)
    
    axes[0].scatter(valid_results['Cp'], valid_results['rho'], alpha=0.5)
    axes[0].set_xlabel('Cp (J/mol·K)')
    axes[0].set_ylabel('Density (mol/m³)')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].scatter(valid_results['RED'], valid_results['Cp'], alpha=0.5)
    axes[1].set_xlabel('RED')
    axes[1].set_ylabel('Cp (J/mol·K)')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].scatter(valid_results['RED'], valid_results['rho'], alpha=0.5)
    axes[2].set_xlabel('RED')
    axes[2].set_ylabel('Density (mol/m³)')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('property_correlations.png', dpi=300, bbox_inches='tight')
    print(f"Saved correlations to 'property_correlations.png'")
    
    # Export to CSV
    valid_results.to_csv('sampled_properties.csv', index=False)
    print(f"\nExported {len(valid_results)} valid results to 'sampled_properties.csv'")
    
    # Show example molecules
    print("\n" + "="*80)
    print("EXAMPLE MOLECULES")
    print("="*80)
    
    for i in range(min(5, len(valid_results))):
        row = valid_results.iloc[i]
        print(f"\nMolecule {i+1} ({row['molecule_type']}):")
        print(f"  Groups: {row['ni_dict']}")
        print(f"  Cp = {row['Cp']:.6f} J/(mol·K)")
        print(f"  rho = {row['rho']:.2f} mol/m³")
        print(f"  RED = {row['RED']:.6f}")
    
    plt.show()

else:
    print("\nNo valid results generated. Try increasing sample size or adjusting constraints.")

print("\n" + "="*80)
print("SAMPLING COMPLETE")
print("="*80)
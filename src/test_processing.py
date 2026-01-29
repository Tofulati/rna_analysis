"""
Test script to verify data processing works correctly
"""
import pickle
import pandas as pd
import numpy as np
from pathlib import Path

def process_gene_table(gene_name, table):
    """
    Process a single gene's table to extract modification data
    """
    gene_stats = {'gene_name': gene_name}
    
    # For each region
    regions = ['UTR_5', 'UTR_3', 'Exon', 'Intron']
    modifications = ['Unmod', 'm6A', 'Inosine']
    
    print(f"\nProcessing {gene_name}:")
    print(f"Table shape: {table.shape}")
    print(f"Features present: {table['Feature'].unique()}")
    
    for region in regions:
        region_key = region.lower().replace('_', '')
        region_data = table[table['Feature'] == region]
        
        if len(region_data) == 0:
            print(f"  {region}: No data found")
            continue
        
        print(f"\n  {region}:")
        
        # Get CPK values
        cpk_1 = region_data['CPK_MR01_1'].values
        cpk_2 = region_data['CPK_MR01_2'].values
        cpk_values = np.concatenate([cpk_1, cpk_2])
        cpm = np.mean(cpk_values[cpk_values > 0]) if len(cpk_values[cpk_values > 0]) > 0 else 0
        gene_stats[f'{region_key}_cpm'] = cpm
        print(f"    CPM: {cpm:.2f}")
        
        # Get modification rates
        for mod_type in modifications:
            mod_data = region_data[region_data['Modification'] == mod_type]
            
            if len(mod_data) > 0:
                mr01_1_str = mod_data['MR01_1'].values[0]
                mr01_2_str = mod_data['MR01_2'].values[0]
                
                try:
                    rate_1 = float(mr01_1_str.split('±')[0].strip())
                    rate_2 = float(mr01_2_str.split('±')[0].strip())
                    rate = (rate_1 + rate_2) / 2
                except:
                    rate = 0
                
                if mod_type == 'Inosine':
                    gene_stats[f'{region_key}_ai_rate'] = rate
                    print(f"    A-to-I rate: {rate:.4f} ({rate*100:.2f}%)")
                elif mod_type == 'm6A':
                    gene_stats[f'{region_key}_m6a_rate'] = rate
                    print(f"    m6A rate: {rate:.4f} ({rate*100:.2f}%)")
    
    # Calculate "either" rates
    for region in regions:
        region_key = region.lower().replace('_', '')
        ai_key = f'{region_key}_ai_rate'
        m6a_key = f'{region_key}_m6a_rate'
        
        if ai_key in gene_stats and m6a_key in gene_stats:
            either = min(gene_stats[ai_key] + gene_stats[m6a_key], 1.0)
            gene_stats[f'{region_key}_either_rate'] = either
    
    # Calculate totals
    rate_keys = ['ai_rate', 'm6a_rate', 'either_rate']
    cpm_values = []
    
    for region in regions:
        region_key = region.lower().replace('_', '')
        cpm_key = f'{region_key}_cpm'
        if cpm_key in gene_stats and gene_stats[cpm_key] > 0:
            cpm_values.append(gene_stats[cpm_key])
        
        for rate_type in rate_keys:
            total_key = f'total_{rate_type}'
            region_rate_key = f'{region_key}_{rate_type}'
            
            if region_rate_key in gene_stats:
                if total_key not in gene_stats:
                    gene_stats[total_key] = []
                gene_stats[total_key].append(gene_stats[region_rate_key])
    
    if cpm_values:
        gene_stats['total_cpm'] = np.mean(cpm_values)
        
    for rate_type in rate_keys:
        total_key = f'total_{rate_type}'
        if total_key in gene_stats and isinstance(gene_stats[total_key], list):
            if len(gene_stats[total_key]) > 0:
                gene_stats[total_key] = np.mean(gene_stats[total_key])
    
    print(f"\n  Total stats:")
    print(f"    Total CPM: {gene_stats.get('total_cpm', 0):.2f}")
    print(f"    Total A-to-I: {gene_stats.get('total_ai_rate', 0)*100:.2f}%")
    print(f"    Total m6A: {gene_stats.get('total_m6a_rate', 0)*100:.2f}%")
    print(f"    Total Either: {gene_stats.get('total_either_rate', 0)*100:.2f}%")
    
    return gene_stats


# Test with AZIN1
database_dir = Path('./database')
test_file = database_dir / 'AZIN1.pkl'

print("="*60)
print("Testing data processing with AZIN1.pkl")
print("="*60)

with open(test_file, 'rb') as f:
    data = pickle.load(f)

print("\nRaw data:")
print(data)

gene_stats = process_gene_table('AZIN1', data)

print("\n" + "="*60)
print("Final gene_stats dictionary:")
print("="*60)
for key, value in sorted(gene_stats.items()):
    if isinstance(value, float):
        print(f"{key:25s}: {value:.4f}")
    else:
        print(f"{key:25s}: {value}")

print("\n" + "="*60)
print("Test complete! If you see reasonable numbers above, the processing is working correctly.")
print("="*60)
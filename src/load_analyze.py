import pickle
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats
import json

class RNAModificationAnalyzer:
    def __init__(self, database_dir='./database'):
        self.database_dir = Path(database_dir)
        self.genes_data = {}
        self.combined_df_mr01_1 = None
        self.combined_df_mr01_2 = None
        self.raw_gene_tables = {}  # Store raw tables for display

    def load_all_genes(self):
        print("loading gene data from pkls")

        pickle_files = list(self.database_dir.glob('*.pkl'))
        print(f"Found {len(pickle_files)} pkls")
        
        # Debug: Check the first file to verify format
        if len(pickle_files) > 0:
            first_file = pickle_files[0]
            print(f"\nDiagnostic check on first file: {first_file.stem}")
            try:
                with open(first_file, 'rb') as f:
                    sample_data = pickle.load(f)
                print(f"  Columns: {sample_data.columns.tolist()}")
                print(f"  Features: {sample_data['Feature'].unique().tolist()}")
                print(f"  Modifications: {sample_data['Modification'].unique().tolist()}")
                print(f"\nFirst few rows:")
                print(sample_data.head())
            except Exception as e:
                print(f"  Error in diagnostic: {e}")

        for pkl in pickle_files:
            gene_name = pkl.stem

            try: 
                with open(pkl, 'rb') as f:
                    data = pickle.load(f)
                self.genes_data[gene_name] = data
                # Store raw table for later export
                self.raw_gene_tables[gene_name] = data.copy()
            except Exception as e:
                print(f"error loading {gene_name}: {e}")

        print(f"\nsucessfully loaded {len(self.genes_data)} genes")
        return self.genes_data
    
    def process_gene_table(self, gene_name, table, sample='MR01_1'):
        """Process gene table for a specific sample (MR01_1 or MR01_2)"""
        gene_stats = {'gene_name': gene_name}

        regions = ['UTR_5', 'UTR_3', 'Exon', 'Intron']
        modifications = ['Unmod', 'm6A', 'Inosine']
        
        # Determine which columns to use based on sample
        if sample == 'MR01_1':
            cpk_col = 'CPK_MR01_1'
            mr_col = 'MR01_1'
            count_col = 'Count_MR01_1'
        else:
            cpk_col = 'CPK_MR01_2'
            mr_col = 'MR01_2'
            count_col = 'Count_MR01_2'

        for region in regions:
            region_key = region.lower().replace('_','')

            region_data = table[table['Feature'] == region]

            if len(region_data) == 0:
                continue

            # Get CPM for this sample only
            cpk_values = region_data[cpk_col].values
            gene_stats[f'{region_key}_cpm'] = np.mean(cpk_values[cpk_values > 0]) if len(cpk_values[cpk_values > 0]) > 0 else 0

            # Initialize rates to 0 for this region
            gene_stats[f'{region_key}_ai_rate'] = 0
            gene_stats[f'{region_key}_m6a_rate'] = 0
            gene_stats[f'{region_key}_unmod_rate'] = 0

            for mod_type in modifications:
                mod_data = region_data[region_data['Modification'] == mod_type]

                if len(mod_data) > 0:
                    mr_str = mod_data[mr_col].values[0]

                    try:
                        # Extract the rate for this sample
                        # Handle both string format "0.084891 ± 0.170510" and numeric
                        if isinstance(mr_str, str):
                            rate = float(mr_str.split('±')[0].strip())
                        else:
                            rate = float(mr_str)
                    except:
                        rate = 0

                    if mod_type == 'Inosine':
                        gene_stats[f'{region_key}_ai_rate'] = rate
                    elif mod_type == 'm6A':
                        gene_stats[f'{region_key}_m6a_rate'] = rate
                    elif mod_type == 'Unmod':
                        gene_stats[f'{region_key}_unmod_rate'] = rate
                    
        # Calculate "either" modification rates for each region
        for region in regions:
            region_key = region.lower().replace('_','')
            ai_key = f'{region_key}_ai_rate'
            m6a_key = f'{region_key}_m6a_rate'

            # Ensure both keys exist
            if ai_key not in gene_stats:
                gene_stats[ai_key] = 0
            if m6a_key not in gene_stats:
                gene_stats[m6a_key] = 0
                
            # Calculate either modification
            gene_stats[f'{region_key}_either_rate'] = min(
                gene_stats[ai_key] + gene_stats[m6a_key], 1.0
            )

        # Calculate total gene metrics (average across regions)
        rate_keys = ['ai_rate', 'm6a_rate', 'either_rate']
        cpm_values = []

        for region in regions:
            region_key = region.lower().replace('_','')
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
        else:
            gene_stats['total_cpm'] = 0

        # Calculate averages for total rates
        for rate_type in rate_keys:
            total_key = f'total_{rate_type}'
            if total_key in gene_stats and isinstance(gene_stats[total_key], list):
                if len(gene_stats[total_key]) > 0:
                    gene_stats[total_key] = np.mean(gene_stats[total_key])
                else:
                    gene_stats[total_key] = 0
            else:
                gene_stats[total_key] = 0

        return gene_stats
    
    def create_combined_dataframe(self):
        """Create separate dataframes for MR01_1 and MR01_2 samples"""
        print("processing gene tables for both samples")

        all_gene_stats_mr01_1 = []
        all_gene_stats_mr01_2 = []

        for gene_name, table in self.genes_data.items():
            try: 
                # Process for MR01_1
                gene_stats_1 = self.process_gene_table(gene_name, table, sample='MR01_1')
                all_gene_stats_mr01_1.append(gene_stats_1)
                
                # Process for MR01_2
                gene_stats_2 = self.process_gene_table(gene_name, table, sample='MR01_2')
                all_gene_stats_mr01_2.append(gene_stats_2)
            except Exception as e:
                print(f"error processing {gene_name}: {e}")

        self.combined_df_mr01_1 = pd.DataFrame(all_gene_stats_mr01_1)
        self.combined_df_mr01_2 = pd.DataFrame(all_gene_stats_mr01_2)
        
        print(f"created MR01_1 dataframe with {len(self.combined_df_mr01_1)} genes")
        print(f"created MR01_2 dataframe with {len(self.combined_df_mr01_2)} genes")
        
        # Debug: Print column names to verify all data was extracted
        print(f"\nMR01_1 columns: {sorted(self.combined_df_mr01_1.columns.tolist())}")
        
        # Debug: Print sample of data for first gene
        if len(self.combined_df_mr01_1) > 0:
            print(f"\nSample gene (MR01_1):")
            first_gene = self.combined_df_mr01_1.iloc[0]
            print(f"  Gene: {first_gene['gene_name']}")
            print(f"  UTR5 - AI: {first_gene.get('utr5_ai_rate', 'N/A'):.4f}, m6A: {first_gene.get('utr5_m6a_rate', 'N/A'):.4f}, Either: {first_gene.get('utr5_either_rate', 'N/A'):.4f}")
            print(f"  UTR3 - AI: {first_gene.get('utr3_ai_rate', 'N/A'):.4f}, m6A: {first_gene.get('utr3_m6a_rate', 'N/A'):.4f}, Either: {first_gene.get('utr3_either_rate', 'N/A'):.4f}")

    def calculate_regression(self, x, y):
        mask = ~(np.isnan(x) | np.isnan(y) | np.isinf(x) | np.isinf(y))
        x_clean = x[mask]
        y_clean = y[mask]

        if len(x_clean) < 2:
            return None
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_clean, y_clean)

        return {
            'slope': slope,
            'intercept': intercept,
            'r_squared': r_value ** 2,
            'p_value': p_value,
            'std_err': std_err,
            'n': len(x_clean)
        }
    
    def create_regression_plot(self, region, mod_type, sample='MR01_1', output_dir='plots'):
        """Create regression plot for a specific sample"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Select the appropriate dataframe
        df_source = self.combined_df_mr01_1 if sample == 'MR01_1' else self.combined_df_mr01_2

        mod_rate_col = f'{region}_{mod_type}_rate'
        cpm_col = f'{region}_cpm'

        if mod_rate_col not in df_source.columns or cpm_col not in df_source.columns:
            print(f"skipping {sample} {region}_{mod_type} - missing columns")
            return None
    
        df = df_source[[mod_rate_col, cpm_col, 'gene_name']].dropna()

        if len(df) < 2:
            print(f"skipping {sample} {region}_{mod_type} - insufficient data")
            return None
        
        x = df[mod_rate_col].values * 100  # Convert to percentage
        y = df[cpm_col].values

        reg_stats = self.calculate_regression(x, y)

        if reg_stats is None:
            return None
        
        # Create figure
        fig, ax = plt.subplots(figsize=(10, 7))

        # Scatter plot
        ax.scatter(x, y, alpha=0.6, s=50, color='#4ECDC4', edgecolors='white', linewidth=0.5)

        # Regression line
        x_line = np.array([x.min(), x.max()])
        y_line = reg_stats['slope'] * x_line + reg_stats['intercept']
        ax.plot(x_line, y_line, 'r--', linewidth=2.5, label=f'R² = {reg_stats["r_squared"]:.4f}')

        # Labels
        mod_labels = {'ai': 'A-to-I', 'm6a': 'm6A', 'either': 'Either Modification'}
        region_labels = {
            'utr5': "5' UTR", 
            'utr3': "3' UTR", 
            'exon': 'Exonic', 
            'intron': 'Intronic', 
            'total': 'Total Gene'
        }
        
        ax.set_xlabel(f'{mod_labels.get(mod_type, mod_type)} (%)', fontsize=13, fontweight='bold')
        ax.set_ylabel('CPM (Counts Per Million)', fontsize=13, fontweight='bold')
        ax.set_title(
            f'{region_labels.get(region, region)} - {mod_labels.get(mod_type, mod_type)} ({sample})',
            fontsize=15,
            fontweight='bold',
            pad=20
        )

        # Statistics box
        stats_text = (
            f"R² = {reg_stats['r_squared']:.4f}\n"
            f"p-value = {reg_stats['p_value']:.2e}\n"
            f"slope = {reg_stats['slope']:.2f}\n"
            f"n = {reg_stats['n']} genes"
        )
        ax.text(
            0.05, 0.95, stats_text,
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'),
            fontsize=11,
            family='monospace'
        )
        
        ax.legend(loc='lower right', fontsize=11)
        ax.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()

        # Save
        filename = f'{sample}_{region}_{mod_type}_regression.png'
        plt.savefig(output_dir / filename, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"created: {filename} (R²={reg_stats['r_squared']:.4f}, n={reg_stats['n']})")
        return reg_stats
    
    def create_all_plots(self, output_dir='plots'):
        """Create all regression plots for both samples"""
        regions = ['utr5', 'utr3', 'exon', 'intron', 'total']
        mod_types = ['ai', 'm6a', 'either']
        samples = ['MR01_1', 'MR01_2']

        all_stats = {}

        print("\n" + "="*60)
        print("GENERATING ALL REGRESSION PLOTS")
        print("="*60)
        
        for sample in samples:
            print(f"\n--- Processing {sample} ---")
            for region in regions:
                for mod_type in mod_types:
                    stats = self.create_regression_plot(region, mod_type, sample, output_dir)
                    if stats:
                        all_stats[f'{sample}_{region}_{mod_type}'] = stats

        print("\n" + "="*60)
        print(f"COMPLETED: Generated {len(all_stats)} plots")
        print("="*60 + "\n")
        
        return all_stats
    
    def export_data_for_react(self, output_dir='output'):
        """Export data in JSON format for React component - both samples with RAW table data"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        
        # Export MR01_1
        genes_mr01_1 = []
        for idx, row in self.combined_df_mr01_1.iterrows():
            gene_name = row['gene_name']
            gene_dict = {
                'id': int(idx + 1),
                'name': gene_name,
                'chromosome': f"chr{np.random.randint(1, 23)}",
            }
            
            # Add aggregated stats
            for col in self.combined_df_mr01_1.columns:
                if col != 'gene_name':
                    value = row[col]
                    if pd.isna(value):
                        value = 0
                    if isinstance(value, (np.integer, np.floating)):
                        value = float(value)
                    gene_dict[col] = value
            
            # Add raw table data if available
            if gene_name in self.raw_gene_tables:
                raw_table = self.raw_gene_tables[gene_name]
                raw_data = []
                
                for _, raw_row in raw_table.iterrows():
                    # Extract the mean value from the MR01_1 column
                    mr_str = raw_row['MR01_1']
                    if isinstance(mr_str, str):
                        try:
                            mr_value = float(mr_str.split('±')[0].strip())
                        except:
                            mr_value = 0
                    else:
                        mr_value = float(mr_str) if not pd.isna(mr_str) else 0
                    
                    raw_data.append({
                        'feature': raw_row['Feature'],
                        'modification': raw_row['Modification'],
                        'count': int(raw_row['Count_MR01_1']) if not pd.isna(raw_row['Count_MR01_1']) else 0,
                        'cpk': float(raw_row['CPK_MR01_1']) if not pd.isna(raw_row['CPK_MR01_1']) else 0,
                        'mr': mr_value
                    })
                
                gene_dict['raw_data'] = raw_data
            
            genes_mr01_1.append(gene_dict)
        
        # Export MR01_2
        genes_mr01_2 = []
        for idx, row in self.combined_df_mr01_2.iterrows():
            gene_name = row['gene_name']
            gene_dict = {
                'id': int(idx + 1),
                'name': gene_name,
                'chromosome': f"chr{np.random.randint(1, 23)}",
            }
            
            # Add aggregated stats
            for col in self.combined_df_mr01_2.columns:
                if col != 'gene_name':
                    value = row[col]
                    if pd.isna(value):
                        value = 0
                    if isinstance(value, (np.integer, np.floating)):
                        value = float(value)
                    gene_dict[col] = value
            
            # Add raw table data if available
            if gene_name in self.raw_gene_tables:
                raw_table = self.raw_gene_tables[gene_name]
                raw_data = []
                
                for _, raw_row in raw_table.iterrows():
                    # Extract the mean value from the MR01_2 column
                    mr_str = raw_row['MR01_2']
                    if isinstance(mr_str, str):
                        try:
                            mr_value = float(mr_str.split('±')[0].strip())
                        except:
                            mr_value = 0
                    else:
                        mr_value = float(mr_str) if not pd.isna(mr_str) else 0
                    
                    raw_data.append({
                        'feature': raw_row['Feature'],
                        'modification': raw_row['Modification'],
                        'count': int(raw_row['Count_MR01_2']) if not pd.isna(raw_row['Count_MR01_2']) else 0,
                        'cpk': float(raw_row['CPK_MR01_2']) if not pd.isna(raw_row['CPK_MR01_2']) else 0,
                        'mr': mr_value
                    })
                
                gene_dict['raw_data'] = raw_data
            
            genes_mr01_2.append(gene_dict)
        
        # Save both as separate JSON files
        json_path_1 = output_dir / 'gene_data_MR01_1.json'
        with open(json_path_1, 'w') as f:
            json.dump(genes_mr01_1, f, indent=2)
        
        json_path_2 = output_dir / 'gene_data_MR01_2.json'
        with open(json_path_2, 'w') as f:
            json.dump(genes_mr01_2, f, indent=2)
        
        # Also save combined structure
        combined_export = {
            'MR01_1': genes_mr01_1,
            'MR01_2': genes_mr01_2
        }
        json_path_combined = output_dir / 'gene_data_combined.json'
        with open(json_path_combined, 'w') as f:
            json.dump(combined_export, f, indent=2)
        
        print(f"Exported {len(genes_mr01_1)} genes (MR01_1) to {json_path_1}")
        print(f"Exported {len(genes_mr01_2)} genes (MR01_2) to {json_path_2}")
        print(f"Exported combined data to {json_path_combined}")
        
        return {
            'MR01_1': genes_mr01_1,
            'MR01_2': genes_mr01_2
        }

    def export_data(self, output_dir='output'):
        """Export data as CSV and JSON for both samples"""
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        # Export CSVs
        self.combined_df_mr01_1.to_csv(output_dir / 'gene_data_MR01_1.csv', index=False)
        self.combined_df_mr01_2.to_csv(output_dir / 'gene_data_MR01_2.csv', index=False)
        
        print(f"Exported MR01_1 gene data to {output_dir / 'gene_data_MR01_1.csv'}")
        print(f"Exported MR01_2 gene data to {output_dir / 'gene_data_MR01_2.csv'}")
        
        # Also export for React
        self.export_data_for_react(output_dir)

    def run_full_analysis(self, output_dir='output'):
        """Run complete analysis pipeline for both samples"""
        output_dir = Path(output_dir)
        
        print("\n" + "="*60)
        print("RNA MODIFICATION ANALYSIS PIPELINE")
        print("="*60 + "\n")
        
        # Load genes
        self.load_all_genes()

        # Process both samples
        self.create_combined_dataframe()

        # Generate all plots
        print("\nGenerating regression plots...")
        stats = self.create_all_plots(output_dir / 'plots')

        # Export data
        print("\nExporting data files...")
        self.export_data(output_dir)

        # Create regression stats summary
        stats_df = pd.DataFrame([
            {
                'sample': k.split('_')[0] + '_' + k.split('_')[1],
                'region': k.split('_')[2],
                'mod_type': k.split('_')[3],
                **v
            }
            for k, v in stats.items()
        ])
        stats_df.to_csv(output_dir / 'regression_stats.csv', index=False)
        
        print(f"\n" + "="*60)
        print("ANALYSIS COMPLETE")
        print("="*60)
        print(f"\nResults saved to: {output_dir.absolute()}/")
        print(f"  - Plots: {(output_dir / 'plots').absolute()}/")
        print(f"  - Data: gene_data_MR01_1.csv, gene_data_MR01_2.csv")
        print(f"  - JSON: gene_data_MR01_1.json, gene_data_MR01_2.json")
        print(f"  - Combined JSON: gene_data_combined.json")
        print(f"  - Stats: regression_stats.csv")
        print(f"\nTotal plots generated: {len(stats)}")
        print("="*60 + "\n")

        return stats
    
if __name__ == '__main__':
    analyzer = RNAModificationAnalyzer(database_dir='./database')

    stats = analyzer.run_full_analysis(output_dir='./output')

    # Print detailed summary
    print("\nREGRESSION STATISTICS SUMMARY:")
    print("-" * 80)
    
    # Group by sample
    for sample in ['MR01_1', 'MR01_2']:
        print(f"\n{sample}:")
        sample_stats = {k: v for k, v in stats.items() if k.startswith(sample)}
        
        for key, value in sorted(sample_stats.items()):
            region = key.split('_')[2]
            mod_type = key.split('_')[3]
            print(f"  {region:8} | {mod_type:6} | R²={value['r_squared']:.4f} | "
                  f"p={value['p_value']:.2e} | n={value['n']:4}")
    
    print("\n" + "="*80)
    print("All files have been generated successfully!")
    print("Upload gene_data_MR01_1.json or gene_data_MR01_2.json to the React app")
    print("="*80)
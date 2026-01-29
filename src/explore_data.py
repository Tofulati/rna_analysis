"""
Explore pickle files to understand the data structure
"""
import pickle
import pandas as pd
from pathlib import Path

# Point to your database directory
database_dir = Path('./database')

# Get list of all pickle files
pickle_files = list(database_dir.glob('*.pkl'))

print(f"Found {len(pickle_files)} pickle files")
print("\nFirst few files:")
for f in pickle_files[:5]:
    print(f"  {f.name}")

# Load one example file to see structure
if pickle_files:
    example_file = pickle_files[0]
    print(f"\n\nExamining: {example_file.name}")
    
    with open(example_file, 'rb') as f:
        data = pickle.load(f)
    
    print(f"\nData type: {type(data)}")
    
    if isinstance(data, pd.DataFrame):
        print(f"\nDataFrame shape: {data.shape}")
        print("\nColumn names:")
        print(data.columns.tolist())
        print("\nFirst few rows:")
        print(data.head())
        print("\nData types:")
        print(data.dtypes)
    else:
        print(f"\nData: {data}")
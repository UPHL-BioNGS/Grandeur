import pandas as pd
import re

def summarize_mash_err(summary_df, mash_err):
    file = mash_err
    print("Adding results for " + file)
    analysis = "mash_err"
    
    with open(file, 'r') as f:
        content = f.read()
    
    pattern = re.compile(
        r"Estimated genome size:\s+([\d.e+]+).*?"
        r"Estimated coverage:\s+([\d.]+).*?"
        r"Writing to\s+(.*?)\.msh", 
        re.DOTALL
    )
    
    matches = pattern.findall(content)
    
    data = []
    for size, coverage, sample_id in matches:
        data.append({
            "sample": sample_id.strip(),
            "genome_size": float(size),
            "coverage": float(coverage)
        })
    
    # Create DataFrame
    new_df = pd.DataFrame(data)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)
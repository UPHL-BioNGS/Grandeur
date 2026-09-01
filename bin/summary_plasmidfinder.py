import json
import pandas as pd

def parse_concatenated_json(json_string):
    """Generator to parse concatenated JSON objects from a string."""
    decoder = json.JSONDecoder()
    json_string = json_string.strip()
    while json_string:
        obj, index = decoder.raw_decode(json_string)
        yield obj
        json_string = json_string[index:].lstrip()

# plasmidfinder : merging relevant rows into one
def summarize_plasmidfinder(summary_df, plasmidfinder):
    file = plasmidfinder
    print("Adding results for " + file)
    analysis = "plasmidfinder"
    with open(file, 'r') as f:
        content = f.read()

    all_data = []
    
    for run in parse_concatenated_json(content):
        # Extract sample name from the JSON value
        execs = run.get("software_executions", {})
        first_exec = next(iter(execs.values())) if execs else {}
        out_json_val = first_exec.get("parameters", {}).get("out_json", "")
        sample_name = out_json_val.split("/")[-2]

        seq_regions = run.get("seq_regions", {})
        
        plasmids = ", ".join([f"{d.get('name', '')} ({d.get('identity', '')})" for d in seq_regions.values()])
        all_data.append({
            "sample": sample_name,
            "plasmid_(identity)": plasmids if plasmids else ""
        })

        
    new_df = pd.DataFrame(all_data)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = new_df.columns.str.replace('Sample', 'sample', regex=True)
    new_df.columns = [x.lower() for x in new_df.columns]
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how = 'left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)
    return summary_df

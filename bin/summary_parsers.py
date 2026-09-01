#!/usr/bin/env python3

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

def parse_file(summary_df, file, delim):
    """Parses standard CSV/TSV metrics and merges them into the summary dataframe."""
    print("Adding results for " + str(file))
    analysis = str(file).split("_")[0]
    new_df = pd.read_csv(file, dtype=str, index_col=False, delimiter=delim)
    new_df = new_df.add_prefix(analysis + "_")
    new_df.columns = new_df.columns.str.replace('Sample', 'sample', regex=True)
    new_df.columns = [x.lower() for x in new_df.columns]
    new_df[analysis + "_sample"] = new_df[analysis + "_sample"].astype(str)
    summary_df = pd.merge(summary_df, new_df, left_on="sample", right_on=analysis + "_sample", how='left')
    summary_df.drop(analysis + "_sample", axis=1, inplace=True)
    return summary_df
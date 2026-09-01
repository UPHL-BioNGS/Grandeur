#!/usr/bin/env python3

import os
import pandas as pd

def create_node(parent=None):
    """Creates a basic dictionary representing a node."""
    return {
        'parent': parent,
        'children': [],
        'name': "",
        'length': 0.0
    }

def parse_newick(newick_str):
    """Parses a newick string into linked dictionaries."""
    newick_str = newick_str.strip().replace('\n', '').replace('\r', '')
    root = create_node()
    current_node = root
    state = 'name'
    buffer = ""

    for char in newick_str:
        if char == '(':
            new_node = create_node(parent=current_node)
            current_node['children'].append(new_node)
            current_node = new_node
            state = 'name'
            buffer = ""
        elif char == ',':
            if state == 'name':
                current_node['name'] = buffer.strip()
            else:
                if buffer.strip():
                    current_node['length'] = float(buffer)
            
            new_node = create_node(parent=current_node['parent'])
            current_node['parent']['children'].append(new_node)
            current_node = new_node
            state = 'name'
            buffer = ""
        elif char == ')':
            if state == 'name':
                current_node['name'] = buffer.strip()
            else:
                if buffer.strip():
                    current_node['length'] = float(buffer)
            
            current_node = current_node['parent']
            state = 'name'
            buffer = ""
        elif char == ':':
            if state == 'name':
                current_node['name'] = buffer.strip()
            state = 'length'
            buffer = ""
        elif char == ';':
            if state == 'name':
                current_node['name'] = buffer.strip()
            elif state == 'length':
                if buffer.strip():
                    current_node['length'] = float(buffer)
            break
        else:
            buffer += char

    return root

def get_tip_distance_stats(newick):
    """Calculates the average, min, and max distances using dictionaries."""
    with open(newick, 'r') as f:
        newick_str = f.read()

    file_basename = os.path.basename(newick)
    root = parse_newick(newick_str)
    
    tips = []
    def find_tips(node):
        if len(node['children']) == 0 and node['name']: 
            tips.append(node)
        for child in node['children']:
            find_tips(child)
            
    find_tips(root)
    
    if len(tips) < 2:
        return pd.DataFrame()
        
    results = {}
    
    for start_tip in tips:
        distances = {}
        visited = set()
        queue = [(start_tip, 0.0)]
        visited.add(id(start_tip))
        
        while queue:
            current, dist = queue.pop(0)
            is_tip = len(current['children']) == 0
            
            if is_tip and id(current) != id(start_tip):
                distances[current['name']] = dist
                
            parent = current['parent']
            if parent and id(parent) not in visited:
                visited.add(id(parent))
                queue.append((parent, dist + current['length']))
                
            for child in current['children']:
                if id(child) not in visited:
                    visited.add(id(child))
                    queue.append((child, dist + child['length']))
                    
        dist_values = list(distances.values())
        results[start_tip['name']] = {
            'average': sum(dist_values) / len(dist_values),
            'min': min(dist_values),
            'max': max(dist_values)
        }
    
    df_data = [{'file': file_basename, 'sample': tip, **metrics} for tip, metrics in results.items()]
    df = pd.DataFrame(df_data)
    df = df.add_prefix(file_basename + "_")
    df['sample'] = df[file_basename + "_sample"]

    return df

def get_snp_distance_stats(filepath):
    """Calculates the average, min, and max distances from a snp-dists matrix."""
    df_matrix = pd.read_csv(filepath, index_col=0)
    file_basename = os.path.basename(filepath)
    results = []
    
    for sample in df_matrix.index:
        distances = df_matrix.loc[sample].drop(str(sample))
        results.append({
            'file': file_basename,
            'sample': sample,
            'average': distances.mean(),
            'min': distances.min(),
            'max': distances.max()
        })
        
    df = pd.DataFrame(results)
    df = df.add_prefix(file_basename + "_")
    df['sample'] = df[file_basename + "_sample"]
    return df
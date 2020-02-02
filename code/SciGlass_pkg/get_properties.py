import pandas as pd
import numpy as np
import requests
import os
import json
from pprint import pprint
import re
os.environ["MAPI_KEY"] = 'gqfAFKShNOzMcXa0'



def get_data(chemsys_formula_id, data_type, prop, MAPI_KEY = 'gqfAFKShNOzMcXa0'):
    import requests, os, json
    import pandas as pd
    os.environ['MAPI_KEY'] = MAPI_KEY
    
    try:
        if isinstance(prop, list) and data_type == 'vasp':
            sub_url = "https://materialsproject.org/rest/v2/materials/%s/%s" %(chemsys_formula_id, data_type)
            request_data = requests.get(sub_url, headers={"X-API-KEY": os.environ["MAPI_KEY"]})
            request_data_df = pd.DataFrame(request_data.json()['response'])
            return request_data_df[prop]
        else:
            if isinstance(prop, list):
                prop = prop.pop()
            sub_url = "https://materialsproject.org/rest/v2/materials/%s/%s/%s" %(chemsys_formula_id, data_type, prop)
            request_data = requests.get(sub_url, headers={"X-API-KEY": os.environ["MAPI_KEY"]})
            return pd.DataFrame(request_data.json()['response'])
    except:
        return pd.DataFrame([[1],[2]])


def simplify_elasticity(property_df):  
    nrow = property_df.shape[0]
    elastic_mat = np.zeros((nrow, 3))
    elastic_df = pd.DataFrame(elastic_mat, columns = ['material_id','K', 'G'])
    elastic_df['material_id'] = property_df['material_id']

    for i in range(nrow):
        elastic_dict = property_df[['elasticity']].iloc[i, 0]
        if elastic_dict is not None:
                elastic_df.iloc[i, 1] = elastic_dict['K_VRH']
                elastic_df.iloc[i, 2] = elastic_dict['G_VRH']
        else:
            elastic_dict = property_df[['elastic_moduli']].iloc[i, 0]
            elastic_df.iloc[i, 1] = elastic_dict['K']
            elastic_df.iloc[i, 2] = elastic_dict['G']
    return elastic_df


def combine_formula(elem_list, necessary_elem, n_elem='all'):
    n = len(elem_list)
    necessary_elem_str = '-'.join(necessary_elem)
    
    import itertools
    if n_elem == 'all':
        all_set = []
        for n_elem_current in range(n):
            all_set.extend(list(itertools.combinations(elem_list, n_elem_current+1)))
    else:
        all_set = itertools.combinations(elem_list, n_elem - 1)
        
    return ['-'.join(elem) + '-' + necessary_elem_str for elem in all_set]




def get_data_combine(elem_list, necessary_elem, data_type, prop, n_elem = 'all', MAPI_KEY = 'gqfAFKShNOzMcXa0'):
    
    combine_list = combine_formula(elem_list = elem_list, necessary_elem = necessary_elem, n_elem = n_elem)
    
    final_df = pd.DataFrame()
    for combine in combine_list:
        temp_df = get_data(combine, data_type = data_type, prop = prop, MAPI_KEY = MAPI_KEY)
        final_df = final_df.append(temp_df, ignore_index=True)    
    return final_df



def get_glass_df(elem_list, necessary_elem, prop, e_above_hull = None, n_elem = 'all', MAPI_KEY = 'gqfAFKShNOzMcXa0'):
    
    glass_df = get_data_combine(elem_list = elem_list, necessary_elem = necessary_elem, data_type = 'vasp', prop = prop)
    if 'elasticity' in prop:
        predicted_elastic_df = get_data_combine(elem_list = elem_list, necessary_elem = necessary_elem, data_type = 'pred', prop = 'elastic_moduli')
        raw_K_G_elastic_df = pd.merge(glass_df[['material_id','elasticity']], predicted_elastic_df, on='material_id')
        
    if e_above_hull is not None:
        glass_df = glass_df[glass_df['e_above_hull']<=e_above_hull]
        
    if  'elasticity' in prop:
        K_G_elastic_df = simplify_elasticity(raw_K_G_elastic_df)
        glass_df =  pd.merge(glass_df.drop(columns=['elasticity']), K_G_elastic_df, on = 'material_id')    
         
    return glass_df

#chemsys_formula_id = 'Na-Mg-O'

def clean_info_df(raw_df):
    raw_df['atom_volume'] = raw_df['volume']/raw_df['nsites']
    raw_df = raw_df.rename(columns = {'formation_energy_per_atom':'formation_energy'})
    raw_df = raw_df.rename(columns = {'pretty_formula':'formula'})
    property_names = ['material_id', 'formula', 'formation_energy', 'band_gap', 'atom_volume' , 'density', 'G', 'K']
    
    clean_df = raw_df[property_names]
    unique_formula_set = set(clean_df['formula'])
    final_df = pd.DataFrame(np.zeros((len(unique_formula_set), len(property_names))), columns = property_names)
    i = 0
    for formula in unique_formula_set:
        temp_sub_df = clean_df[clean_df['formula'] == formula]
        final_df.iloc[i,:] = temp_sub_df[temp_sub_df['formation_energy'] == np.min(temp_sub_df['formation_energy'])].values[0]  
        i += 1
    return final_df






prop = ['material_id', 'full_formula', 'pretty_formula', 'formation_energy_per_atom',
        'e_above_hull', 'band_gap','density', 'volume', 'nsites', 'elasticity']
elem_list = ['Ca', 'Mg', 'Al', 'Si', 'Na', 'K', 'B']
necessary_elem = ['O']

a = get_glass_df(elem_list, necessary_elem, prop, e_above_hull = 0.005, n_elem = 'all', MAPI_KEY = 'gqfAFKShNOzMcXa0')
b = clean_info_df(a)
b.to_csv('7_info.csv')



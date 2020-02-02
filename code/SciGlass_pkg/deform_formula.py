from get_properties import *


def unit_atom(atoms_list, element_dict, times = 1):
    if len(atoms_list) == 1:
        return element_dict
    else:
        element_dict[atoms_list[0]] = 1 * times
        return unit_atom(atoms_list[1:], element_dict, times)

def raw_count_atom(formula, element_dict, times = 1):
    atoms = re.findall(r'[A-Za-z]+', formula)
    num_atoms = re.findall(r'\d+', formula)
    num_atoms = [int(num) for num in num_atoms]

    for atom in atoms:
        temp = re.findall(r'[A-Z][a-z]*', atom)
        if len(temp) >= 2:
            element_dict = unit_atom(temp, element_dict, times)
            atoms[atoms.index(atom)] = temp[-1]
        else: 
            continue
    
    for i in range(len(num_atoms)):
        element_dict[atoms.pop(0)] += num_atoms.pop(0) * times
    if atoms:
        element_dict[atoms[0]] += 1 * times
    
    return element_dict

def count_atom(formula, element_dict):
    sub_formulas = re.findall(r"\((.*?)\)", formula)
    sub_num = re.findall(r"\(.*?\)(?:(\d+)?)", formula)
    sub_num = [int(num) if num.isdigit() else 1 for num in sub_num]
    for sub_f in sub_formulas:
        element_dict = raw_count_atom(sub_f, element_dict, times = sub_num.pop(0))

    sub_formulas = re.split(r" ?\([^)]+\)", formula)
    sub_formulas = [sub_f for sub_f in sub_formulas if len(sub_f) != 0 and not sub_f.isdigit()]
    for sub_f in sub_formulas:
        element_dict = raw_count_atom(sub_f, element_dict)
    return element_dict



prop = ['material_id', 'full_formula', 'pretty_formula', 'formation_energy_per_atom',
        'e_above_hull', 'band_gap','density', 'volume', 'nsites', 'elasticity']
#elem_list = ['Ca', 'Mg', 'Al', 'Si', 'Na', 'K', 'B']
elem_list = ['Ca', 'Mg', 'Al', 'Si', 'Na', 'K']
necessary_elem = ['O']

a = get_glass_df(elem_list, necessary_elem, prop, e_above_hull = 0.01, n_elem = 'all', MAPI_KEY = 'gqfAFKShNOzMcXa0')
info_df = clean_info_df(a)

#oxide_names = ['CaO','MgO','Al2O3','SiO2', 'Na2O', 'K2O', 'B2O3']
oxide_names = ['CaO','MgO','Al2O3','SiO2', 'Na2O', 'K2O']
compound_mat = np.zeros((info_df.shape[0], len(oxide_names)))
element_mat = np.zeros((info_df.shape[0], len(oxide_names)+1))
i = 0
for formula in info_df['formula']:
    #element_dict = {'Ca':0, 'Mg':0, 'Al':0, 'Si':0, 'Na':0 ,'K': 0, 'B':0, 'O':0}
    element_dict = {'Ca':0, 'Mg':0, 'Al':0, 'Si':0, 'Na':0 ,'K': 0, 'O':0}
    element_dict = count_atom(formula, element_dict)  
    #compound_mat[i,:] = np.array(list(element_dict.values())[:len(oxide_names)])/np.array([1,1,2,1,2,2,2])
    compound_mat[i,:] = np.array(list(element_dict.values())[:len(oxide_names)])/np.array([1,1,2,1,2,2])
    element_mat[i,:] = np.array(list(element_dict.values()))[:(len(oxide_names)+1)]
    i += 1
    
compound_df = pd.DataFrame(compound_mat, columns = oxide_names)
element_df = pd.DataFrame(element_mat, columns = element_dict.keys())
#print(element_df.head(10))
#element_df.to_csv('Deform.csv', header = True, index = None)
#bool_idx = np.equal(np.sum(element_df[['Ca', 'Mg', 'Al', 'Si', 'Na', 'K', 'B']].values*np.array([[1,1,1.5,2,0.5,0.5,1.5]]), axis = 1), element_df[['O']].values.squeeze())
bool_idx = np.equal(np.sum(element_df[['Ca', 'Mg', 'Al', 'Si', 'Na', 'K']].values*np.array([[1,1,1.5,2,0.5,0.5]]), axis = 1), element_df[['O']].values.squeeze())
b = pd.concat([info_df, compound_df,element_df], axis = 1)
b[bool_idx].to_csv('info.csv', header = True, index = None)



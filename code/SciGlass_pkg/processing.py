import numpy as np
import pandas as pd

def normalize(x_df, axis = 1):
	'''
	Normalize a pd.DataFrame by row and column.

	---------------
	Parameters:
		x_df: pd.DataFrame
			A pd.DataFrame ready to normalized. 
		axis: integer
			0 or 1, if 0, normalize df by column, or by row. (Default 1)

	---------------
	Return: 
		A normalized pd.DataFrame.
	'''	

	x_mat = x_df.values
	sum_up = np.sum(x_mat, axis = axis)
	normlized_mat = x_mat/sum_up.reshape((-1,1))
	return pd.DataFrame(normlized_mat, index = x_df.index, columns = x_df.columns)



def sphere_trans(x_df):
	if "DataFrame" in str(type(x_df)):
		x = x_df.values
	else:
		raise TypeError('Input should be pandas.DataFrame.')
	
	nrow, ncol = x.shape
	y = np.zeros((nrow, ncol-1))
	
	for col_idx in np.arange(ncol-1):
		y[:, col_idx] = np.arctan2(x[:, col_idx+1:].sum(1)**.5, x[:,col_idx]**.5)
	
	return pd.DataFrame(y, index = x_df.index)


def property_name_mod(info_df):
	'''
	Unify the names of properties.

	---
	Parameters:
		info_df: pd.DataFrame
			A pd.DataFrame containing properties
	---
	Return:
		A pd.DataFrame with new properties' names
	'''
	property_name = ['Formation Energy (eV)', 'Band Gap (eV)','Atomic volume', 'Density (gm/cc)', 'G (GPa)', 'K (GPa)']
	new_property_name = ["formation_energy", "band_gap", "atom_volume", "density" , "G", "K"]

	name_pair = dict(zip(property_name, new_property_name))
	
	return info_df.rename(columns = name_pair)





def create_dist_df(data_df, info_df, formula_name = 'Formula'):
	'''
	Compute the a pd.DataFrame of distances.

	---------------
	Parameters:
		data_df: pd.DataFrame 
			The raw dataset. (only included composition)
		info_df: pd.DataFrame 
			The dataset of solid compounds.
		formula_name: string 
			The column name of the formula column.
	---------------
	Return:
		A pd.DataFrame of Euclidean Distances between data points 
		and solid compounds. (shape: #data * #compounds)
	'''

	if formula_name not in info_df.columns:
		formula_name = 'formula'

	#if not all(['O' in oxide for oxide in data_df.columns]):
		#raise ValueError('data_df is not valid.')

	n_data = data_df.shape[0]  # Number of data points
	n_solid = info_df.shape[0] # Number of solid compounds
	n_base = data_df.shape[1]  # Number of basic oxides (e.g. CaO, SiO2, MgO, ...)
	
	dist_mat = np.zeros([n_data, n_solid]) # Initialize the distance-matrice
	for i in range(n_data):
		delta = data_df.iloc[i,:]-info_df[data_df.columns]
		dist_mat[i,:] = np.linalg.norm(delta , ord=2, axis=1)

	return pd.DataFrame(dist_mat, index = data_df.index, columns = info_df[formula_name])



def create_weights(energy, dist, width, intercept, top = 21, normalized = True):
	'''
	Create a bunch of weighting factors for ONLY ONE data point.
	
	---------------
	Parameters:
		energy: pd.Series 
			Energy engaging in constructing weighting factors. (optional energy: Formation Energy or Cohesive Energy)
		dist: pd.Series 
			A dataframe with one row of distances.
		width: float 
			The bandwidth engaging in constructing weighting factors. 
		top: integer
			Only choose [top] solid compound with first [top] largest weighting factor.
		normalized: bool 
			Whether normalize the weighting factors or not. (Default True)
	
	---------------
	Return:
		A pd.Series and the part of exponent.
	'''     
	
	# we do not use some tricks to make all kinds of cases have the same frame
	# In the Class NewDescriptors, they will be tidied up???
	d = dist.values[np.isfinite(dist).values]
	h = np.abs(energy.values[np.isfinite(dist).values])

	# Transforming distance and energy into those within [0, 1]
	norm_d = (d-np.min(d)) / (np.max(d)-np.min(d))
	norm_h = (h - np.min(h)) / (np.max(h)-np.min(h))
	fraction = (norm_d)/(norm_h + intercept)

	raw_weight = np.exp(-(fraction/width)**2)

	# Choose Top [top] solid compound: set other weighting factors to zero
	#raw_weight = np.where((-raw_weight).argsort().argsort() < top, raw_weight, 0)

	# Normalized weighting factors
	if normalized:
		return raw_weight/np.sum(raw_weight), fraction
	else: 
		return raw_weight, fraction
			


def create_weight_df(data_df, dist_df, info_df, width, intercept, formula_name = 'Formula', 
					 energy_name = 'formation_energy', top = 21, normalized = True):
	'''
	Create a weighting factors dataframe for all the data points.
	
	---------------
	Parameters:
		data_df: pd.DataFrame 
			A raw dataset (only included composition).
		info_df: pd.DataFrame 
			A dataset of solid compounds.
		dist_df: pd.DataFrame 
			A distance-dataframe produced by create_dist_df function.
		formula_name: string 
			The column name of the formula column.
		energy_name: string 
			The column name of the energy column.
		top: integer 
			Only choose [top] solid compound with first [top] largest weighting factor.
		normalized: bool 
			Whether normalize the weighting factors or not. (Default True)
	
	---------------
	Return:
		A DataFrame of weighting factor and the DataFrame of exponents. (shape: #data * #solids)
	'''
	
	if formula_name not in info_df.columns:
		formula_name = 'formula'

	#if energy_name not in info_df.columns:
	#    energy_name = 'formation_energy'

	weight_mat = np.zeros((dist_df.shape[0], info_df.shape[0]))
	power_mat = np.zeros((dist_df.shape[0], info_df.shape[0]))
	
	for i in range(data_df.shape[0]):
		#idx = np.isfinite(dist_df.iloc[i,:]).values
		weight_mat[i, :], power_mat[i, :] = create_weights(info_df[energy_name], dist_df.iloc[i,:],
														 normalized = normalized, intercept = intercept,
														 width = width, top = top)
		
	return pd.DataFrame(weight_mat, index = data_df.index, columns = info_df[formula_name]), pd.DataFrame(power_mat, index = data_df.index, columns = info_df[formula_name])





import numpy as np
import pandas as pd
from .processing import *

class PropertyDescriptor:
	'''
	Construct descriptors associated with properties
	Four types of stats: weighted mean / weighted standard deviation/ 
						 weighted absoluted deviation/ maximum with largest weight
	'''

	def __init__(self, property_df, weight_df):
		'''
		Initialize the PropertyDescriptor object
		Parameters:
			property_df: pd.DataFrame - a dataframe containing properties of solid compounds
			weight_df: 
		'''
		if 'DataFrame' not in str(type(property_df)):
		    raise TypeError('Properties is not a data frame.')
		if 'DataFrame' not in str(type(weight_df)):
		    raise TypeError('Properties is not a data frame.')
		if property_df.shape[0] != weight_df.shape[1]:
			raise ValueError('Property_df does NOT match with weight_df')

		self.property = property_name_mod(property_df)
		self.weight = weight_df

		self.n_property = self.property.shape[1]
		self.n_data = self.weight.shape[0]

		self.long_property = np.tile(self.property.values.T, (self.n_data, 1))
		self.long_weight = np.repeat(self.weight.values, self.n_property, axis = 0)

	def mean(self):
	    property_mat = self.property.values.reshape((self.property.values.shape[0], -1))
	    mean_mat = np.matmul(self.weight.values, property_mat)
	    col_names = ['avg_' + p for p in self.property.columns]
	    return pd.DataFrame(mean_mat, index = self.weight.index, columns = col_names)

	def max(self):
	    max_idx = np.argmax(self.weight.values, axis = 1)
	    max_mat = self.property.values[max_idx,:]
	    col_names = ['max_' + p for p in self.property.columns]
	    return pd.DataFrame(max_mat, index = self.weight.index, columns = col_names)


	def sd(self):
	    property_mat = self.property.values.reshape((self.property.values.shape[0], -1))
	    mean_mat = np.matmul(self.weight.values, property_mat)

	    n_property = self.n_property
	    n_data = self.n_data

	    long_property = self.long_property
	    long_weight = self.long_weight
	    long_mean = mean_mat.flatten('C').reshape((-1,1))

	    col_names = ['sd_' + p for p in self.property.columns]
	    sd_long_vec = np.sqrt(np.sum((long_property - long_mean)**2 * long_weight, axis = 1))
	    sd_mat = sd_long_vec.reshape((-1, n_property))
	    return pd.DataFrame(sd_mat, index = self.weight.index, columns = col_names)


	def ad(self):
	    property_mat = self.property.values.reshape((self.property.values.shape[0], -1))
	    mean_mat = np.matmul(self.weight.values, property_mat)

	    n_property = self.n_property
	    n_data = self.n_data

	    long_property = self.long_property
	    long_weight = self.long_weight
	    long_mean = mean_mat.flatten('C').reshape((-1,1))

	    col_names = ['ad_' + p for p in self.property.columns]
	    ad_long_vec = np.sum(np.abs(long_property - long_mean) * long_weight, axis = 1)
	    ad_mat = ad_long_vec.reshape((-1, n_property))
	    return pd.DataFrame(ad_mat, index = self.weight.index, columns = col_names)




class LiquidDescriptor:

	def __init__(self, data_df):
		self.data = data_df

	def phasedisorder(self):
		log_x = np.log(self.data.values, where = (self.data.values!=0))
		phase_vec = np.sum(self.data.values * log_x ,axis = 1)
		return pd.DataFrame(phase_vec, index = self.data.index, columns = ['Phase Disorder'])

	def sd(self):
		def unit_sd(x):
			return np.std(x[np.nonzero(x)[0]])
		sd_vec = self.data.apply(lambda x: unit_sd(x), axis = 1)
		return pd.DataFrame(sd_vec, index = self.data.index, columns = ['sd_liquid'])

	def norm(self, order):
		order_df = pd.DataFrame()
		n_len = len(order)
		for i in range(n_len):
			name = 'L'+str(order[i])
			order_df[name] = np.linalg.norm(self.data.values, ord = order[i], axis = 1)
		order_df.index = self.data.index
		return order_df


class NewDescriptors:

	def __init__(self, data_file, info_file, transform = False):
	    
		self.raw_data_df = pd.read_csv(data_file)
		oxide_name = [col_name for col_name in self.raw_data_df.columns if 'O' in col_name]

		self.raw_data_df = pd.DataFrame(normalize(self.raw_data_df[oxide_name]),
		                            index = self.raw_data_df.index, columns = oxide_name)
		
		self.data_df = self.raw_data_df

		self.info_df = property_name_mod(pd.read_csv(info_file))

		if transform:
			self.data_df = sphere_trans(self.data_df)
			new_base_df = sphere_trans(self.info_df[oxide_name])
			self.info_df[new_base_df.columns] = new_base_df


		self.formula_name = 'Formula'

		if 'Formula' not in self.info_df.columns:
			self.formula_name = 'formula'

#		if 'formula' not in self.info_df.columns:
#			raise ValueError("'Formula' or 'formula' column does not in the info data frame")

		self.info_df[oxide_name] = normalize(self.info_df[oxide_name])
		self.dist_df = create_dist_df(self.data_df, self.info_df)

	def weight_df(self, width, intercept, formula_name = 'Formula', energy_name = 'formation_energy', normalized = True, top = None):
		if not top:
		    top = self.info_df.shape[0]
		weight, _ = create_weight_df(self.data_df, self.dist_df, self.info_df, width = width, intercept = intercept,formula_name = formula_name, 
		              energy_name = energy_name, top = top, normalized = normalized)
		return weight

	def property_df(self, property_name, width, intercept, **args):
		weight = self.weight_df(width, intercept, **args)
		#idx = [np.where(self.info_df[self.formula_name] == solid)[0][0] for solid in weight.columns]
		#sub_info_df = self.info_df.iloc[:, idx]
		property_descriptors = PropertyDescriptor(property_df = self.info_df[property_name], weight_df = weight)
		df = [property_descriptors.mean(), property_descriptors.sd(), 
		      property_descriptors.ad(), property_descriptors.max()]
		return pd.concat(df, axis = 1)

	def enthropy_df(self):
		liquid_df = LiquidDescriptor(data_df = self.raw_data_df)
		return liquid_df.phasedisorder()

	def all_df(self, property_name, **args):
	    df = [self.enthropy_df(), self.property_df(property_name = property_name, **args)]
	    return pd.concat(df, axis = 1)







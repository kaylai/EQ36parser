import pandas
import csv
import re
import numpy
import matplotlib.pyplot as plt
from scipy.interpolate import spline

###SAVE SOME CONSTANTS FOR LATER USE###
#Molecular weights dict
MW_dict = {'MAGNETIT':111.69, 'BRUCITE':58.325, 'CLINOCHL':555.8245, 'PYROPE':403.13}

#Molecular weights of solid solutions
def calcMWolivine(Fa):
	MW_olivine = 28.0855+4*(16)+2*(55.845*Fa+24.305*(1-Fa))
	return MW_olivine

def calcMWopx(Fe):
	MW_opx = 2*(55.845*Fe + 24.305*(1-Fe)) + 2*28.0855 + 6*16
	return MW_opx

def calcMWcpx(He):
	MW_cpx = 40.078+(55.845*He + 24.305*(1-He)) + 28.0855*2 + 16*6
	return MW_cpx

#initialize rows
rows = []

#create column names that we actually want as single strings
my_col_names = ['log-zi', 'time', 'log-days']
bad_col_names = ['log zi', 'time, d', 'log days']
more_bad_col_names = ['log  zi', 'time,    d', 'log  days']

with open('tab.csv', 'rb') as file:
	file.seek(0)
	reader = csv.reader(file)
	for line in file:
		rows.append(line)

	startingline_minerals = [i for i in range(len(rows)) if str('log of moles of product') in rows[i]]
	endingline_minerals = [i for i in range(len(rows)) if str('log of destroyed moles of reactants') in rows[i]]
	startingline_minerals[0] += 5
	endingline_minerals[0] -=3

	startingline_solidsolutions = [i for i in range(len(rows)) if str('solid solution product compositions') in rows[i]]
	endingline_solidsolutions = [i for i in range(len(rows)) if str('log of moles of product minerals') in rows[i]]
	startingline_solidsolutions[0] += 2
	endingline_solidsolutions[0] -= 4


mineral_rows = rows[startingline_minerals[0]:endingline_minerals[0]]
solidsolution_rows = rows[startingline_solidsolutions[0]:endingline_solidsolutions[0]]

mineral_rows_clean_headers = []
solidsolution_rows_clean_headers = []

for row in mineral_rows:
	for i in range(len(my_col_names)): #replaces col names with spaces to names without spaces
		row = row.replace(bad_col_names[i], my_col_names[i]) 
		row = row.replace(more_bad_col_names[i], my_col_names[i])
	row = row.lstrip(' ') #removes leading whitespace on each row
	row = re.sub(' +',' ',row) #turns multiple spaces within each row into a single space
	mineral_rows_clean_headers.append(row)

for row in solidsolution_rows:
	row = row.replace('log zi', 'log-zi')
	row = row.lstrip(' ')
	row = re.sub(' +', ' ',row)
	solidsolution_rows_clean_headers.append(row)

df = pandas.DataFrame(mineral_rows_clean_headers)
c = df.columns[0]
df = df[c].str.split(' ', expand=True) #split rows into columns, where columns separated by a single space
df = df.replace('\n', numpy.nan, regex=True) #replaces \n values with NaN's (easier to scrub later)
df = df.replace('None', numpy.nan, regex=True) #replaces 'None' values with NaN's (easier to scrub later)
df = df.dropna(axis=1, how='all')

df.columns = df.iloc[0]
df = df.reindex(df.index.drop(0))

df2 = pandas.DataFrame(solidsolution_rows_clean_headers)
c2 = df2.columns[0]
df2 = df2[c2].str.split(' ', expand=True)
df2 = df2.replace('\n', numpy.nan, regex=True)
df2 = df2.replace('None', numpy.nan, regex=True)
df2 = df2.dropna(axis=1, how='all')

df2.columns = df2.iloc[0]
df2 = df2.reindex(df2.index.drop(0))

#Get header rows
header_list = [idx for idx, row in df.iterrows() if row['log-zi'] == str('log-zi')] #gets all rows that are headers of a new table
header_list.append(0) #add first chunk of values
header_list.sort() #sorts list numerically

header_list_ss = [idx for idx, row in df2.iterrows() if row['log-zi'] == str('log-zi')]
header_list_ss.append(0)
header_list_ss.sort()

#make a dict from a list
header_dict = dict((index, pandas.DataFrame()) for index in header_list)
header_dict_ss = dict((index, pandas.DataFrame()) for index in header_list_ss)

iter_range = len(header_list) - 1
iter_range_ss = len(header_list_ss) - 1

final_list_value = header_list[-1]
final_dataframe_index = df.index[-1]
final_chunk_len = final_dataframe_index - final_list_value + 1

final_list_value_ss = header_list_ss[-1]
final_dataframe_index_ss = df2.index[-1]
final_chunk_len_ss = final_dataframe_index_ss - final_list_value_ss + 1

#create a dict where key = rownumber for header rows and value = pandas dataframe of all data corresponding to that header
for key, value in header_dict.iteritems():
	for i in range(iter_range):
		if header_list[i] == key:
			header_dict[key] = df[header_list[i]:header_list[i+1]-1]
		if key == final_list_value:
			header_dict[final_list_value] = df[-final_chunk_len+1:]

for key, value in header_dict.iteritems():
	if key == 0:
		pass
	else:
		header_dict[key].columns = df.iloc[key-1]
		header_dict[key] = header_dict[key].iloc[1:]
		#header_dict[key] = header_dict[key].reindex(header_dict[key].index.drop(key))
	header_dict[key] = header_dict[key].dropna(axis=1, how='all')
	column_names = header_dict[key].columns.get_values().tolist()
	header_dict[key] = header_dict[key].dropna(subset=[column_names[0]])

for key, value in header_dict_ss.iteritems():
	for i in range(iter_range_ss):
		if header_list_ss[i] == key:
			header_dict_ss[key] = df2[header_list_ss[i]:header_list_ss[i+1]-1]
		if key == final_list_value_ss:
			header_dict_ss[final_list_value_ss] = df2[-final_chunk_len_ss+1:]

for key, value in header_dict_ss.iteritems():
	if key == 0:
		pass
	else:
		header_dict_ss[key].columns = df2.iloc[key-1]
		header_dict_ss[key] = header_dict_ss[key].iloc[1:]
	header_dict_ss[key] = header_dict_ss[key].dropna(axis=1, how='all')
	column_names_ss = header_dict_ss[key].columns.get_values().tolist()

#combine separate dataframes into one and sort by index value
cleaned_df = pandas.concat((header_dict[key] for key, value in header_dict.iteritems()), sort=False)
cleaned_df = cleaned_df.sort_index(axis=0, ascending=True)

cleaned_df2 = pandas.concat((header_dict_ss[key] for key, value in header_dict_ss.iteritems()), sort=False)
cleaned_df2 = cleaned_df2.sort_index(axis=0, ascending=True)

#remove any rows with no numeric values
cleaned_df = cleaned_df[pandas.to_numeric(cleaned_df['log-zi'], errors='coerce').notnull()]
cleaned_df2 = cleaned_df2[pandas.to_numeric(cleaned_df2['log-zi'], errors='coerce').notnull()]

#save numeric values as float instead of text
for col in cleaned_df.columns[0:]:
	cleaned_df[col] = cleaned_df[col].astype(float)
for col in cleaned_df2.columns[0:]:
	cleaned_df2[col] = cleaned_df2[col].astype(float)

#create a dataframe with values as moles of minerals (instead of log moles)
minmoles_df = pandas.DataFrame()

for col in cleaned_df.columns[3:]:
	minmoles_df[col+'_moles'] = 10**(cleaned_df[col])

minmoles_df = minmoles_df.fillna(0)
minmoles_df['moles total minerals'] = sum((minmoles_df[col] for col in minmoles_df.columns))
minmoles_df = minmoles_df.drop_duplicates()

#create a dataframe with values of molar proportions of each mineral
molpropmin_df = pandas.DataFrame()

for col in minmoles_df.columns:
	if col == 'moles total minerals':
		pass
	else:
		molpropmin_df[col.replace('_moles', '_prop')] = minmoles_df[col] / minmoles_df['moles total minerals']

molpropmin_df['log-zi'] = cleaned_df['log-zi']

#save numeric values as float instead of text for excel output
for col in molpropmin_df.columns:
	molpropmin_df[col] = molpropmin_df[col].astype(float)

#set indices to log-zi
cleaned_df = cleaned_df.set_index('log-zi')
cleaned_df2 = cleaned_df2.set_index('log-zi')
molpropmin_df = molpropmin_df.set_index('log-zi')

#Remove log-zi=-999 value row
cleaned_df = cleaned_df.iloc[1:]
#molpropmin_df = molpropmin_df.iloc[1:]


##------Remove duplicate indices------##
cleaned_df = cleaned_df.drop_duplicates(subset=None, keep='first', inplace=False) #delete duplicate total row
cleaned_df2 = cleaned_df2.drop_duplicates(subset=None, keep='first', inplace=False) #delete duplicate total row
minmoles_df = minmoles_df.drop_duplicates(subset=None, keep='first', inplace=False) #delete duplicate total row
molpropmin_df = molpropmin_df.drop_duplicates(subset=None, keep='first', inplace=False) #delete duplicate total row

#Remove duplicate indices with distinct row vals
new_index = cleaned_df.index.values.tolist()
new_index2 = cleaned_df2.index.values.tolist()

for i in range(1,len(new_index)):
	if new_index[i] == new_index[i-1]:
		new_index[i] += 0.0001

cleaned_df['new_index'] = new_index
cleaned_df = cleaned_df.set_index('new_index')
cleaned_df.index.name = 'log-zi'

molpropmin_df['new_index'] = new_index #molpropmin has identical index to df
molpropmin_df = molpropmin_df.set_index('new_index') 
molpropmin_df.index.name = 'log-zi'

minmoles_df['new_index'] = new_index #minmoles has identical index to df
minmoles_df = minmoles_df.set_index('new_index')
minmoles_df.index.name = 'log-zi'

for i in range(1,len(new_index2)):
	if new_index2[i] == new_index2[i-1]:
		new_index2[i] += 0.0001

cleaned_df2['new_index2'] = new_index2
cleaned_df2 = cleaned_df2.set_index('new_index2')
cleaned_df2.index.name = 'log-zi'
#--------------------------------------#


##------Calculate grams of minerals------##
#Create a dataframe combining minmoles and solid solution values (df2)
moles_and_ss_df = pandas.concat([minmoles_df, cleaned_df2], axis=1, sort=False, join='outer')	
moles_and_ss_df = moles_and_ss_df.rename(columns= lambda x: x.replace('_moles', ''))
moles_and_ss_df = moles_and_ss_df.fillna(0)
moles_and_ss_df = moles_and_ss_df.astype(float)

#create a dataframe with grams of minerals
mingrams_df = pandas.DataFrame()

#create a list of all solid solution phases present
solid_solutions_dict = {}
if 'FAYALITE' in list(moles_and_ss_df.columns.values):
	solid_solutions_dict["OLIVINE"] = ["FAYALITE", calcMWolivine]
if 'FERROSILIT' in list(moles_and_ss_df.columns.values):
	solid_solutions_dict["ORTHOPYR"] = ["FERROSILIT", calcMWopx]
if 'HEDENBERGI' in list(moles_and_ss_df.columns.values):
	solid_solutions_dict["CLINOPYR"] = ["HEDENBERGI", calcMWcpx]

for col in moles_and_ss_df.columns.values:
	if col in MW_dict:
		mingrams_df[col] = moles_and_ss_df[col] * MW_dict[col]
	if col in solid_solutions_dict:
		mingrams_df[col] = moles_and_ss_df[col] * solid_solutions_dict[col][1](moles_and_ss_df[solid_solutions_dict[col][0]])

mingrams_df["Total grams"] = sum((mingrams_df[col] for col in mingrams_df.columns))
mingrams_df["Total kg"] = mingrams_df["Total grams"] / 1000.0
#----------------------------------------#

##------Calculate Fe3+/Fetot and fl/rock ratio and save to new df------##
#create new df
redox_output_df = pandas.DataFrame()

#set log-zi as the index
redox_output_df['new_index'] = new_index
redox_output_df = redox_output_df.set_index('new_index')
redox_output_df.index.name = 'log-zi'

#calculate and save some values
redox_output_df["moles_Fe3+"] = minmoles_df["MAGNETIT_moles"]
redox_output_df["moles_Fe2+"] = (
								minmoles_df["MAGNETIT_moles"] + 
								2*minmoles_df["OLIVINE_moles"] * moles_and_ss_df["FAYALITE"] + 
								minmoles_df["ORTHOPYR_moles"] * moles_and_ss_df["FERROSILIT"] + 
								minmoles_df["CLINOPYR_moles"] * moles_and_ss_df["HEDENBERGI"]
								)
redox_output_df["Fe3+/Fetot_molar"] = redox_output_df["moles_Fe3+"] / (redox_output_df["moles_Fe3+"] + redox_output_df["moles_Fe2+"])
redox_output_df["fl/rk wt ratio"] = 1 / mingrams_df["Total kg"]

#---------------------------------------------------------------------#


##------Write all data to an excel file------##
writer = pandas.ExcelWriter('cleaned_tab.xlsx', engine='xlsxwriter')
cleaned_df.to_excel(writer, sheet_name='logMoles Minerals')
cleaned_df2.to_excel(writer, sheet_name='Solid Solutions')
molpropmin_df.to_excel(writer, sheet_name='Mineral proportions')
writer.save()
#--------------------------------------------#


###------PLOTTING------####
cleaned_df = cleaned_df.drop(['time', 'log-days'], axis=1)

###BELOW: some trys toward producing smoothed curves###
# index_values = numpy.array(cleaned_df.index)
# index_values = index_values.astype(numpy.float)
# for item in range(len(index_values)):
# 	if index_values[item] == 0.0:
# 		index_values[item] = 0.00000001

# print index_values
# print index_values.min()
# print index_values.max()

# xnew = numpy.linspace(index_values.min(),index_values.max(),300) #300 represents number of points to make between min and max index values
# column_names = cleaned_df.columns.values

# smoothed_df = pandas.DataFrame(index=xnew)
# for column_name in column_names:
# 	smoothed_df[column_name] = spline(index_values,cleaned_df[column_name],xnew) #TODO getting divide by zero error

cleaned_df.plot()
plt.savefig('minerals.png')

molpropmin_df.plot()
plt.savefig('modes.png')

fig, ax1 = plt.subplots()
redox_output_df.plot(ax=ax1, x='fl/rk wt ratio', y='Fe3+/Fetot_molar', logx=True)
ax1.invert_xaxis()
ax1.axhspan(0.20, 0.30, alpha=0.7, color='#fff0aa')
plt.savefig('Fe23.png')
#-------------------------#

###------PRINTING OUTPUT TO TERMINAL------###
flrk_arc_values = {}
mineralogy_arc_values = {}
mineralogy_by_mineral = {}
for index, row in redox_output_df.iterrows():
	if row["Fe3+/Fetot_molar"] < 0.30 and row["Fe3+/Fetot_molar"] > 0.20:
		flrk_arc_values[index] = row["fl/rk wt ratio"]

header_list = list(molpropmin_df.columns.get_values())
for index, row in molpropmin_df.iterrows():
	if index in flrk_arc_values.keys():
		for i in range(len(header_list)):
			mineralogy_by_mineral[header_list[i]] = row[header_list[i]]
		mineralogy_arc_values[index] = mineralogy_by_mineral

print "\n"
print "RESTULS OF THE MODEL RUN:"
print "When Fe3+/Fetot is between 0.2-0.3 in the solids:"
if len(flrk_arc_values) > 1:
	print "     fl/rk ratio = " + str(round(min(flrk_arc_values.values()), 2)) + "-" + str(round(max(flrk_arc_values.values()), 2))
elif len(flrk_arc_values) == 1:
	print "     fl/rk ratio = " + str(round(flrk_arc_values.values()[0], 2))

print "\n"
print "     Mineralogy: "
for key, value in mineralogy_arc_values.iteritems():
	for key, value in value.iteritems():
		if value > 0.0:
			print "     " + str(key) + " = " + str(round(value, 2))

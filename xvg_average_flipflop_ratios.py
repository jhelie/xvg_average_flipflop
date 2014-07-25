#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_average_flipflop_ratios', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
***************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_average_flipflop
***************************************************

[ DESCRIPTION ]
 
This script calculate the average ratios profile for basic and hydrophobic contacts
made by several flipflopping lipids.

It calculates the avg and (unbiased) std dev and can deal with NaN.

The ratios are calculated for each supplied file, then averaged with an equal weights
(each flipflop is considered as informative as the other).

NB:
the script may give out a warning 'return np.mean(x,axis)/factor', it's ok. it's just
scipy warning us that there were only nans on a row, the result will be a nan as we
expect (see this thread: https://github.com/scipy/scipy/issues/2898).

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o	flipflop_avg	: name of outptut file
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs='+', dest='xvgfilenames', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["flpflop_avg"], help=argparse.SUPPRESS)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if len(args.xvgfilenames) == 1:
	print "Error: only 1 data file specified."
	sys.exit(1)
	
for f in args.xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
		
	global nb_rows
	global nb_cols
	global distance
	global data_ratios_basic
	global data_ratios_hphob

	nb_rows = 0
	nb_cols = 0
		
	for f_index in range(0,len(args.xvgfilenames)):
		#display progress
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		
		#get file content
		filename = args.xvgfilenames[f_index]
		with open(filename) as f:
			lines = f.readlines()
		
		#determine legends and nb of lines to skip
		tmp_nb_rows_to_skip = 0
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
		
		#get data
		tmp_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(tmp_data)[0]
			distance = np.zeros((nb_rows, 1))													#distance from cluster
			data_ctct_basic = np.zeros((nb_rows, len(args.xvgfilenames)))
			data_ctct_hphob = np.zeros((nb_rows, len(args.xvgfilenames)))
			data_ctct_total = np.zeros((nb_rows, len(args.xvgfilenames)))
		else:
			if np.shape(tmp_data)[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
		#check that each file has the same number of columns
		if f_index == 0:
			nb_cols = np.shape(tmp_data)[1]
		else:
			if np.shape(tmp_data)[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[1]) + " data columns, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
		#check that each file has the same first column
		if f_index == 0:
			distance[:,0] = tmp_data[:,0]
		else:
			if not np.array_equal(tmp_data[:,0],distance[:,0]):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#store data
		data_ctct_basic[:,f_index] = tmp_data[:,1]
		data_ctct_hphob[:,f_index] = tmp_data[:,2]
		data_ctct_total[:,f_index] = tmp_data[:,3]

		#remove non meaningful zeros (i.e. not sampled)
		data_ctct_total[data_ctct_total == 0] = np.nan
		data_ctct_basic[np.isnan(data_ctct_total)] = np.nan
		data_ctct_hphob[np.isnan(data_ctct_total)] = np.nan

		#calculate ratios
		tmp_div = np.zeros((nb_rows, 1))
		tmp_div[:,0] = np.nansum(data_ctct_basic + data_ctct_hphob, axis = 1)
		tmp_div[np.isnan(tmp_div)] = 1
		data_ratios_basic = data_ctct_basic / tmp_div
		data_ratios_hphob = data_ctct_basic / tmp_div
		
	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_avg():													#DONE

	global data_ratios_basic_avg
	global data_ratios_basic_std
	global data_ratios_hphob_avg
	global data_ratios_hphob_std
	
	data_ratios_basic_avg = scipy.stats.nanmean(data_ratios_basic, axis = 1)
	data_ratios_basic_std = scipy.stats.nanstd(data_ratios_basic, axis = 1)
	
	data_ratios_hphob_avg = scipy.stats.nanmean(data_ratios_hphob, axis = 1)
	data_ratios_hphob_std = scipy.stats.nanstd(data_ratios_hphob, axis = 1)

	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_average_flipflop_ratios v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	output_xvg.write("# -> weight = " + str(len(args.xvgfilenames)) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label \"distance from cluster z axis (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"average contact profile by type during flipflop\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 4\n")
	output_xvg.write("@ s0 legend \"basic contacts (avg)\"\n")
	output_xvg.write("@ s1 legend \"basic contacts (std)\"\n")
	output_xvg.write("@ s2 legend \"hydrophobic contacts (avg)\"\n")
	output_xvg.write("@ s3 legend \"hydrophobic contacts (std)\"\n")
	
	#data
	for r in range(0, nb_rows):
		results = str(distance[r,0])
		results += "	" + "{:.6e}".format(data_ratios_basic_avg[r,0]) + "	" + "{:.6e}".format(data_ratios_basic_std[r,0]) + "	" + "{:.6e}".format(data_ratios_hphob_avg[r,0]) + "	" + "{:.6e}".format(data_ratios_hphob_std[r,0])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\n\nWriting average file..."
calculate_avg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)

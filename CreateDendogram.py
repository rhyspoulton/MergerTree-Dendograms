#This imports the python3 print function if this code is being run in python2
from __future__ import print_function

from ETFio import LoadETFCatalogue
from createPlotArrays import createPlotArrays
from ReadConfig import  plotOptions
from plotDendogram import plotDendogram
import sys
import numpy as np
import os

if(len(sys.argv)<5):
	raise SystemExit("Incorrect number of arguments parsed.\n \tUsage: CreateDendogram.py <ETF file> <Num plot> <output directory> plot_config.cfg\n")


# SelIndex = 150597  # 139360

#VEL new bad: 29337
#VEL new good: 150598
#VEL paramset-2 good: 139336
#VEL paramset-2 bad: 33883
#VEL old good: 184263
#VEL old bad: 21890
#AHF good: 52
#AHF bad: 42
#Rock good: 57
#Rock bad 41


try:
	nPlot = int(sys.argv[2])
except ValueError:
	raise SystemExit("Please parse a int for the number to be plotted")

#Make the output directory if it does to exist
outdir = sys.argv[3]
if(os.path.isdir(outdir)==False):
	os.mkdir(outdir)

#Get the plot options
plotOpt = plotOptions(sys.argv[4],outdir)

#Load the data from the ETF catalogue
opt,halodata = LoadETFCatalogue(sys.argv[1],plotOpt)

#Get the indexes of the nPlot largest branches
endSnap = "Snap_%03d" %(opt.endSnap)
indexes=np.argsort(halodata[endSnap]["Mass"])[::-1][:nPlot][::-1]

#Loop over all the indexes producing dendograms
for SelIndex in indexes:

	plotData,branchIndicator,depthIndicator,sortIndx,mainBranchIDs = createPlotArrays(opt,plotOpt,halodata,SelIndex)

	#Check if there is anythin to plot
	if(len(branchIndicator)==0):
		continue

	plotDendogram(plotOpt,plotData,depthIndicator,branchIndicator,sortIndx,halodata["Snap_%03d" %opt.endSnap]["HaloID"][SelIndex],mainBranchIDs)

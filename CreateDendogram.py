#This imports the python3 print function if this code is being run in python2
from __future__ import print_function

import sys
sys.path.append(sys.argv[0].replace("CreateDendogram.py",""))
from ETFio import LoadETFCatalogue
from createPlotArrays import createPlotArrays
from ReadConfig import  plotOptions
from plotDendogram import plotDendogram
import numpy as np
import os

if(len(sys.argv)<5):
	raise SystemExit("Incorrect number of arguments parsed.\n \tUsage: CreateDendogram.py <ETF file> <Num plot> <output directory> plot_config.cfg\n")


#Make the output directory if it does to exist
outdir = sys.argv[3]
if(os.path.isdir(outdir)==False):
	os.mkdir(outdir)

#Get the plot options
plotOpt = plotOptions(sys.argv[4],outdir)

#Update the number of dendograms to be plotted
try:
	plotOpt.nPlot = int(sys.argv[2])
except ValueError:
	raise SystemExit("Please parse a int for the number to be plotted")

#Load the data from the ETF catalogue and indexes which the dendograms are to be plotted
opt,halodata,indexes = LoadETFCatalogue(sys.argv[1],plotOpt)

#Loop over all the indexes producing dendograms
for SelIndex in indexes:

	plotData,branchIndicator,depthIndicator,sortIndx,mainBranchIDs = createPlotArrays(opt,plotOpt,halodata,SelIndex)

	#Check if there is anythin to plot
	if(len(branchIndicator)==0):
		continue

	plotDendogram(plotOpt,plotData,depthIndicator,branchIndicator,sortIndx,halodata["Snap_%03d" %opt.endSnap]["HaloID"][SelIndex],mainBranchIDs)
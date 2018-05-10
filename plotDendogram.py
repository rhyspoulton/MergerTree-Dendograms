import matplotlib.pyplot as plt 
from matplotlib import gridspec
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection 
import numpy as np
import time




def setColData(colData,mainBranchIDs=[]):
	"""
	This function processes the colour data into colours for plotting 

	"""

	# MainBranchIDs = np.fromfile()

	tmpColData = np.empty(colData.shape,dtype="str")

	# tmpColData[colData==0] = "b"
	# tmpColData[colData>0] = "r"

	haloSel = colData==-1
	tmpColData[colData==-1] = "b"

	mainSubSel = np.isin(colData,mainBranchIDs)

	tmpColData[mainSubSel] = "r"

	tmpColData[np.invert(mainSubSel | haloSel)]="g"

	return tmpColData

def setColDataColourbar(plotOpt,fig,colData,plotWidth,plotHeight,plotMargins):
	"""
	Function process the colour data for plotting and also create a colour bar
	
	"""

	colData[colData>0] = colData[colData>0]

	norm = plt.Normalize(vmin = np.min(colData), vmax = np.max(colData))
	colors = plt.cm.ScalarMappable(norm= norm,cmap = plotOpt.colMap)

	tmpColData = colors.to_rgba(colData)

	colors.set_array(np.linspace(np.min(colData),np.max(colData),100))

	cbaxes = fig.add_axes([1 -0.9*plotMargins[1]/plotWidth , plotMargins[2]/plotHeight, 0.3 * plotMargins[1]/plotWidth, 1 - plotMargins[3]/plotHeight - plotMargins[2]/plotHeight]) 

	cb = fig.colorbar(colors, cax=cbaxes)
	cb.ax.tick_params(labelsize=20) 

	cb.set_label(plotOpt.cbarLabel,fontsize = 20)

	return tmpColData






def setsizeData(plotOpt,sizeData):
	"""
	This function sets the scale for the size data

	"""
	if(plotOpt.logged):
		minSize = np.min(sizeData[sizeData>0])

		sizeData = sizeData/minSize

		sizeData[sizeData>0] = np.log10(sizeData[sizeData>0])

		sizeData=np.multiply(plotOpt.sizePoint,sizeData,casting="unsafe",dtype="float64")


	else:
		#Find maximum for the data 
		maxSize = np.max(sizeData)

		#Normalize the size data to this max size
		sizeData = sizeData/maxSize

		#Multiply by the size of point desired
		sizeData=np.multiply(plotOpt.sizePoint,sizeData,casting="unsafe",dtype="float64")


	return sizeData



def plotDendogram(plotOpt,plotData,depthIndicator,branchIndicator,sortIndx,SelID,mainBranchIDs):

	print("Plotting the tree")
	start = time.time()

	#Remove branches that live for less than the minimum number of snaps and remove branches deeper than maxdepth
	if(plotOpt.plotSubhaloBranches):
		sel = np.where(((np.count_nonzero(plotData["xposData"],axis=0) >= plotOpt.minNsnapsExist) & (depthIndicator<plotOpt.maxdepth)) | (depthIndicator<2))[0]
	else:
		sel = np.where(((np.count_nonzero(plotData["xposData"],axis=0) >= plotOpt.minNsnapsExist) & (depthIndicator<plotOpt.maxdepth) & (depthIndicator>-1)) | (depthIndicator<2))[0]
	plotData["xposData"] = plotData["xposData"][:,sel]
	plotData["SizeData"] = plotData["SizeData"][:,sel]
	plotData["ColData"] = plotData["ColData"][:,sel]
	if(plotOpt.overplotdata):
		plotData["OverPlotData"] = plotData["OverPlotData"][:,sel]
	depthIndicator = depthIndicator[sel]
	branchIndicator = branchIndicator[sel]
	sortIndx = sortIndx[sel]
	
	#Lets remove branches that never appear in the plot
	sel  = [False if(all(plotData["xposData"][:,i][plotData["xposData"][:,i]>0]>=plotOpt.plotNumRvir)) else True for i in range(plotData["xposData"].shape[1]) ] 
	plotData["xposData"] = plotData["xposData"][:,sel]
	plotData["SizeData"] = plotData["SizeData"][:,sel]
	plotData["ColData"] = plotData["ColData"][:,sel]
	if(plotOpt.overplotdata):
		plotData["OverPlotData"] = plotData["OverPlotData"][:,sel]
	depthIndicator = depthIndicator[sel]
	branchIndicator = branchIndicator[sel]
	sortIndx = sortIndx[sel]



	#Now we will see if any of the deeper branches merging branch has been removed, if so we will remove those from the plot
	sel = (np.in1d(branchIndicator,sortIndx)) | (branchIndicator<0)
	plotData["xposData"] = plotData["xposData"][:,sel]
	plotData["SizeData"] = plotData["SizeData"][:,sel]
	plotData["ColData"] = plotData["ColData"][:,sel]
	if(plotOpt.overplotdata):
		plotData["OverPlotData"] = plotData["OverPlotData"][:,sel]
	depthIndicator = depthIndicator[sel]
	branchIndicator = branchIndicator[sel]
	sortIndx = sortIndx[sel]

	plotData["xposData"][plotData["xposData"]==-1.0] = plotOpt.plotNumRvir 

	maxBranchSize = np.amax(plotData["SizeData"],axis=0)

	#Find the number of branches left and the first and final snapshot in the data
	numBranches = plotData["xposData"].shape[1] 
	if(numBranches==1):
		raise Exception("After requiring that the halos exist for %i snaps there is only the main branch left, \nThis can also occur because the units are not correct please check this." %plotOpt.minNsnapsExist)

	fsnap = plotData["xposData"].shape[0] -1 +plotOpt.snapoffset
	isnap = np.where(np.count_nonzero(plotData["xposData"],axis=1) > 0)[0][0] -10 + plotOpt.snapoffset

	#If the number of branches is greater than the plotOpt.on set it to the plotOpt.on
	if(numBranches>plotOpt.maxNumBranches):
		numBranches = plotOpt.maxNumBranches

	#Set the height and width of the plot in inches
	plotHeight = 16
	plotWidth = 1.5*numBranches+1

	#Intialize the plot and the subplots within it
	fig = plt.figure(figsize=(plotWidth,plotHeight))
	gs = gridspec.GridSpec(1,numBranches,width_ratios=[plotOpt.numSubplotsMain if(i==0) else 1 for i in range(numBranches)])

	axes = [fig.add_subplot(gs[i]) for i in range(numBranches)]

	branchSel = np.where((depthIndicator<=1) & (branchIndicator<0))[0][:4]
	
	if(plotOpt.insetPlot):
		if(numBranches>=4):
			
			insetData = plotData["SizeData"][:,branchSel]
			insetCol =["green","orange","purple","pink"]
			insetLS = ["-","--","-.",":"]
		else:
			print("The inset plot cannot be plotted as there are less than 4 branches")

	plotData["SizeData"] = setsizeData(plotOpt,plotData["SizeData"])

	#Set the plot margins of the plot in inches [left,right,bottom,top] depending if there is a colorbar or not on the plot
	if(plotOpt.plotColorBar):
		plotMargins = [1.3,2.5,1.7,1.0]

	else:
		plotMargins = [1.3,0.5,1.7,1.0]


	# Set the col and size data from the functions
	if(plotOpt.plotColorBar):
		plotData["ColData"]=setColDataColourbar(plotOpt,fig,plotData["ColData"],plotWidth,plotHeight,plotMargins)
	else:
		plotData["ColData"] = setColData(plotData["ColData"],mainBranchIDs)

	norm = plt.Normalize(vmin = 1, vmax = 5)
	prog_colors = plt.cm.ScalarMappable(norm= norm,cmap = "Purples_r")

	norm = plt.Normalize(vmin = -1, vmax =5)
	sub_colors = plt.cm.ScalarMappable(norm= norm,cmap = "Purples_r")
	maxDist = np.max(plotData["xposData"][:,0]) + np.max(plotData["SizeData"][:,0])

	ibranchSel=0


	for i in range(numBranches):

		#Find where the data is non-zero for this branch these will be the snapshots
		snaps=np.where(plotData["xposData"][:,i]>0)[0]

		axes[i].set_ylim(isnap,fsnap)
		yticks=np.arange((isnap+9) //10 *10,fsnap+10,10)
		axes[i].set_yticks(yticks)
		axes[i].grid(True,axis="y")
		axes[i].tick_params(axis='both', which='major', labelsize=20)


		#Add subplot and put the data on it
		# ax = fig.add_subplot(gs[i])

		x = plotData["xposData"][:,i][snaps]
		y = snaps +plotOpt.snapoffset

		if(plotOpt.marker=="line"):

			patches = []

			for snap in snaps:


				if(plotOpt.logged):
					width = plotData["SizeData"][:,i][snap] + np.log10(maxDist/plotOpt.plotNumRvir) if(i>0) else plotData["SizeData"][:,i][snap]
				else:
					width = plotData["SizeData"][:,i][snap] * maxDist/plotOpt.plotNumRvir if(i>0) else plotData["SizeData"][:,i][snap]
				xi = plotData["xposData"][:,i][snap] - width/2
				yi = snap - 0.5 +plotOpt.snapoffset
				height = 1
				patch = Rectangle([xi,yi],width,height)
				patches.append(patch)

			lines = PatchCollection(patches,linewidth=0,color=plotData["ColData"][:,i][snaps])
			# lines.set_array(plotData["ColData"][:,i][snaps].decode("UTF-8"))

			axes[i].add_collection(lines)

			
		elif(plotOpt.marker=="circle"):
			axes[i].scatter(x,y,s=plotData["SizeData"][:,i][snaps],c=plotData["ColData"][:,i][snaps],alpha = 0.5,marker ="o" )
		else:
			raise SystemExit("Invalid marker style, please select from either 'line' or 'circle'")
		#Set the y data and the ticks
		
		if(plotOpt.overplotdata):
			for snap in snaps:
				axes[i].text(0,snap - 0.5 + plotOpt.snapoffset,plotOpt.overplotFormat %(plotData["OverPlotData"][:,i][snap]),fontsize=8)
			

		if((i!=(numBranches-1)) & (i>0)):
			if(branchIndicator[i]==0):
				axes[i].scatter(plotOpt.plotNumRvir/2.0,isnap+5,s=2000,color=prog_colors.to_rgba(depthIndicator[i]))
			elif((depthIndicator[i]==1) & (depthIndicator[i+1]>1)):
				axes[i].fill_between([0.5,plotOpt.plotNumRvir],isnap,isnap+10,color=prog_colors.to_rgba(depthIndicator[i]))
			elif((depthIndicator[i]>1) & (depthIndicator[i+1]<2)):
				axes[i].fill_between([0.0,plotOpt.plotNumRvir-0.5],isnap,isnap+10,color=prog_colors.to_rgba(depthIndicator[i]))
			elif((depthIndicator[i]>1) & (depthIndicator[i+1]>1)):
				axes[i].fill_between([0.0,plotOpt.plotNumRvir],isnap,isnap+10,color = prog_colors.to_rgba(depthIndicator[i]))
			elif((depthIndicator[i]==-1) & (depthIndicator[i+1]>1)):
				axes[i].fill_between([0.5,plotOpt.plotNumRvir],isnap,isnap+10,color="yellow")
			elif(depthIndicator[i]==-1):
				axes[i].fill_between([0.5,plotOpt.plotNumRvir-0.5],isnap,isnap+10,color="yellow")
			else:
				axes[i].fill_between([0.5,plotOpt.plotNumRvir-0.5],isnap,isnap+10,color=prog_colors.to_rgba(depthIndicator[i]))
		elif(i>0):
			if(depthIndicator[i]>1):
				axes[i].fill_between([0.0,plotOpt.plotNumRvir-0.5],isnap,isnap+10,color = prog_colors.to_rgba(depthIndicator[i]))
			elif(depthIndicator[i]<0):
				axes[i].fill_between([0.5,plotOpt.plotNumRvir-0.5],isnap,isnap+10,color="yellow")
			else:
				axes[i].fill_between([0.5,plotOpt.plotNumRvir-0.5],isnap,isnap+10,color=prog_colors.to_rgba(depthIndicator[i]))




		#If in the middle of the plot add the y label
		if(i==np.rint(numBranches/2)):
			axes[i].set_xlabel(plotOpt.xLabel,fontsize=30)
			axes[i].text(0,fsnap+7,plotOpt.sizeLabel,fontsize=30)

		#Add a dashed line ar 1 Rvir if not on the main branch and set the ylim to 2.5
		if(i>0):
				axes[i].text(plotOpt.plotNumRvir/2.0 - 0.1,isnap-20,"%i" %i,fontsize=30, weight = 'bold')
				axes[i].axvline(1,ls="dashed")
				axes[i].set_xlim(0,plotOpt.plotNumRvir)
				axes[i].set_xticks(np.arange(1,plotOpt.plotNumRvir,1))
				axes[i].set_yticklabels(['']*len(yticks))
		else:

				axes[i].set_ylabel("Snapshot",fontsize=30)
				axes[i].set_xlabel("Euclidean distance [Mpc]",fontsize=30)
				axes[i].margins(x=0.0)
				# axes[i].set_xlim(0,maxDist)
				for label in axes[i].xaxis.get_ticklabels()[::2]: label.set_visible(False)

		#Connect up the the points
		for j in range(len(x)-1):
				if((numBranches>=4) & (i in branchSel) & (plotOpt.insetPlot)):
					axes[i].plot(x,y,lw=4,ls=insetLS[ibranchSel] ,color=insetCol[ibranchSel])
				else: 
					axes[i].plot((x[j],x[j+1]),(y[j],y[j+1]), color="black")
		
		if(branchIndicator[i]==-3):
			x = x[-1]
			y = y[-1]
			c = plotData["ColData"][:,i][snaps][-1]
			s = 500
			axes[i].scatter(x,y,c=c,s=s,marker="x",linewidths=3)
		elif(branchIndicator[i]==-4):
			x = x[-1]
			y = y[-1]
			c = plotData["ColData"][:,i][snaps][-1]
			s = 500
			axes[i].scatter(x,y,c=c,s=s,marker="^",linewidths=3)
	


		axes[i].tick_params(axis='both', which='major', labelsize=24)
		ax2 = axes[i].twiny()
		# ax3 = ax.twiny()
		ax2.set_xlim(axes[i].get_xlim())
		# ax3.set_xlim(ax.get_xlim())
		ax2.set_xticks([])
		# ax3.set_xticks([])
		ax2.set_xlabel(plotOpt.maxSizeFormat %(maxBranchSize[i]),fontsize = 22)
		# new_fixed_axis = ax3.get_grid_helper().new_fixed_axis
		# ax3.axis["top"] = new_fixed_axis(loc="top",axes=ax3,offset=(0, offset))
		# ax3.axis["top"].toggle(all=True)
		# ax3.axis["top"].label.set_fontsize(7)
		# ax3.set_xlabel("%i" %RootTails[i])

		if((numBranches>=4) & (i in branchSel)):
			ibranchSel+=1

	if((numBranches>=4) & (plotOpt.insetPlot)):
		if(plotOpt.plotColorBar):
			inset = fig.add_axes([0.694,0.215,0.26,0.2])
		else:
			inset = fig.add_axes([0.694,0.215,0.29,0.2])

		for sel in range(4):
			snap = np.where(insetData[:,sel]>0)[0] 
			inset.semilogy(snap+plotOpt.snapoffset,insetData[:,sel][snap],lw=4,ls=insetLS[sel],color = insetCol[sel])


		inset.set_xlabel("Snapshot",fontsize=30)
		inset.set_ylabel(plotOpt.sizeLabel,fontsize=30)
		inset.tick_params(axis='both', which='major', labelsize=24)

	#Adjust the margins of the plot
	plt.subplots_adjust(wspace=0,left=plotMargins[0]/plotWidth,right=1-plotMargins[1]/plotWidth,bottom=plotMargins[2]/plotHeight,top = 1 - plotMargins[3]/plotHeight )


	#Save the file
	outfile = plotOpt.outdir+"/%i_mergertree." %(SelID)+plotOpt.fileDesc+".png"
	plt.savefig(outfile)

	#Close the plot
	plt.close("all")

	print("Done plotting the tree in",time.time() - start)










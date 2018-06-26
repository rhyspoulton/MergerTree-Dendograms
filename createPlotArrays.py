
####################################################################################

# Functions to get all RootProgenitors and create the plotting arrays

# By: Rhys Poulton

####################################################################################



import numpy as np 
from collections import deque
import cosFuncs
import time
from scipy.interpolate import interp1d

def getRootProgenitors(opt,treedata,SelIndex):

	start = time.time()

	#Set the snapname for the final snapshot in the simulation
	endSnapName = "Snap_%03d" %opt.endSnap

	#Set the main branches Root Descendant and Root Progeintor
	MainBranchRootDescedant = treedata[endSnapName]["RootDescendant"][SelIndex]
	MainBranchRootProgenitor = treedata[endSnapName]["RootProgenitor"][SelIndex]

	#Get the final and initial snapshots for this main branch
	MainBranchisnap = int(MainBranchRootProgenitor/opt.HALOIDVAL)
	MainBranchfsnap = int(MainBranchRootDescedant/opt.HALOIDVAL)
	MainBranchNsnaps = MainBranchfsnap - MainBranchisnap

	if(MainBranchisnap<opt.startSnap):
		raise SystemExit("The main branches root progenitor is at snapshot %i but starting snapshot for the plotting is at %i so the full branch cannot be plotted. Please adjust the startSnap for plotting" %(MainBranchisnap,opt.startSnap))

	if(MainBranchNsnaps<2):
		print("The branch with Root Descendant %d only exist for a single snapshot so cannot be plotted." %MainBranchRootDescedant)
		return [],[]

	#A container which has all the RootProgenitors of interest for this Tree
	TreeRootProgenitors = deque([MainBranchRootProgenitor])

	print("Finding tree RootProgenitors")
	#Move down the tree to see which halos point to the MainBranchRootDescedant for this tree
	for snap in range(opt.startSnap,opt.endSnap):

		snapname = "Snap_%03d" %(snap)

		#Find which halo at this snapshot point to the MainBranchRootDescedant
		sel = np.where(treedata[snapname]["RootDescendant"]==MainBranchRootDescedant)[0]

		if(len(sel)>0):
			#Add it to the list of the Tree Root Progenitors
			TreeRootProgenitors.extend(treedata[snapname]["RootProgenitor"][sel].astype(np.int64))

	# Make these unique so they only appear once
	TreeRootProgenitors = np.unique(TreeRootProgenitors)

	#find the original root tail of interest and place it at zero
	sel=np.where(TreeRootProgenitors==MainBranchRootProgenitor)[0]
	if (sel!=0):
		TreeRootProgenitors[0],TreeRootProgenitors[sel] = TreeRootProgenitors[sel],TreeRootProgenitors[0]

	#This will be the identifier if a branch is a main or subbranch(-1) sub-subbranch (>-1) or anything deeper (-3)
	ProgenBranchIndicator = -1 * np.ones(len(TreeRootProgenitors),dtype=np.int32)

	# # Now we have the unique tree root progenitors we can now move up the tree to find the depth of the tree
	# for ibranch,RootProgenitor in enumerate(TreeRootProgenitors[1:]):

	# 	#First extract the haloID of the rootprogenitor, its index and snapshot
	# 	haloID = RootProgenitor
	# 	haloIndex = int(haloID%opt.HALOIDVAL - 1)
	# 	haloSnap = int(haloID/opt.HALOIDVAL)

	# 	snapname = "Snap_%03d" %haloSnap

	# 	#Extract its descendant, index and snapshot
	# 	descID = treedata[snapname]["Descendant"][haloIndex]
	# 	descIndex = int(descID%opt.HALOIDVAL - 1)
	# 	descSnap = int(descID/opt.HALOIDVAL)

	# 	#Parameter to keep track of what depth we are currently at
	# 	depth = 1

	# 	while(True):

	# 		snapname = "Snap_%03d" %descSnap
	# 		#Extract the denscendants progenitor
	# 		progenID = treedata[snapname]["Progenitor"][descIndex]

	# 		#See if the descedants progenitor points back down to the descendant's haloID
	# 		if(progenID!=haloID):

	# 			#If not then extract its Root Progenitor of the branch it merges with
	# 			BranchRootProgen = treedata[snapname]["RootProgenitor"][descIndex]

	# 			#If this doesn't point to the main branch, keep track of which one it points to 
	# 			if(BranchRootProgen!=MainBranchRootProgenitor):
	# 				sel = np.where(TreeRootProgenitors==BranchRootProgen)[0]
	# 				ProgenBranchIndicator[ibranch+1] = sel

	# 			break

	# 		#Otherwise continue to move up the branch extracting its ID, index and snapshot
	# 		haloID = descID
	# 		descID = treedata[snapname]["Descendant"][descIndex]
	# 		descIndex = int(descID%opt.HALOIDVAL-1)
	# 		descSnap =  int(descID/opt.HALOIDVAL)


	#Now we need to identify all subhalos across the lifetime of the main Branch plus include its host if the main branch is itself a subhalo
	SubhaloBranchRootProgenitors = deque()

	snapname = endSnapName
	haloID = MainBranchRootDescedant

	haloIndex = int(haloID%opt.HALOIDVAL -1)
	haloSnap = int(haloID/opt.HALOIDVAL)

	#Store the progenitor main branch so we can identify all subhalos
	progenID = treedata[snapname]["Progenitor"][haloIndex]
	progenIndex = int(progenID%opt.HALOIDVAL -1)
	progenSnap = int(progenID/opt.HALOIDVAL)

	#Store the host root progen if the main branch is ever a subhalo
	MainBranchHostRootProgenitorIDs = []

	print("Found",len(TreeRootProgenitors),"Tree RootProgenitors")

	while(True):

		#Lets identify all halos which have the main branch and store their Root Progenitors so we can check they are also a progenitor so have been plotted twice
		sel = np.where((treedata[snapname]["HostHaloID"]==haloID) & (treedata[snapname]["RootDescendant"]!=MainBranchRootDescedant))

		if(len(sel)>0):
			SubhaloBranchRootProgenitors.extend(treedata[snapname]["RootProgenitor"][sel].astype(np.int64))
		

		if(treedata[snapname]["HostHaloID"][haloIndex]!=-1):
			HostHaloID = treedata[snapname]["HostHaloID"][haloIndex]
			hostIndex = int(HostHaloID%opt.HALOIDVAL-1)
			SubhaloBranchRootProgenitors.extend([treedata[snapname]["RootProgenitor"][hostIndex].astype(np.int64)])
			MainBranchHostRootProgenitorIDs.extend([treedata[snapname]["RootProgenitor"][hostIndex].astype(np.int64)])
		if(haloID == progenID): break

		#Now lets move to the progenitor
		haloID = progenID
		haloIndex = progenIndex
		haloSnap = progenSnap
		snapname = "Snap_%03d" %haloSnap

		if(haloSnap<opt.startSnap):
			raise SystemExit("Going to a halo at snapshot %i but the start snapshot for plotting is %i so it cannot be plotted, please adjust the startSnap" %(haloSnap,opt.startSnap))

		#Find its progenitor now
		progenID = treedata[snapname]["Progenitor"][haloIndex]
		progenIndex = int(progenID%opt.HALOIDVAL - 1)
		progenSnap = int(progenID/opt.HALOIDVAL)

	print("Found",len(SubhaloBranchRootProgenitors),"Subhalo RootProgenitors")

	
	#Only if found any halos which were subhalos of the main branch
	if(len(SubhaloBranchRootProgenitors)>0):

		#Lets make all these subhalo root progenitors unique
		SubhaloBranchRootProgenitors = np.unique(SubhaloBranchRootProgenitors)

		#Set the indicator for the subhalo branches -2 
		SubhaloBranchesIndicator = -2 * np.ones(len(SubhaloBranchRootProgenitors),dtype=np.int32)

		#Now concatenate both the TreeRootProgenitors and the SubhaloBranchRootProgenitors, plus their Indicators
		AllRootProgenitors = np.concatenate([TreeRootProgenitors,SubhaloBranchRootProgenitors[::-1]])
		AllBranchIndicators = np.concatenate([ProgenBranchIndicator,SubhaloBranchesIndicator[::-1]])

		#Make them all unique so a branch does not appear more than once
		_,uniqueIndexes = np.unique(AllRootProgenitors,return_index=True)

		uniqueIndexes=np.sort(uniqueIndexes)

		AllRootProgenitors = AllRootProgenitors[uniqueIndexes]
		AllBranchIndicators = AllBranchIndicators[uniqueIndexes]

		#find the original root tail of interest and place it at zero
		sel=np.where(AllRootProgenitors==MainBranchRootProgenitor)[0]
		if (sel!=0):
			AllRootProgenitors[0],AllRootProgenitors[sel] = AllRootProgenitors[sel],AllRootProgenitors[0]


		if(len(MainBranchHostRootProgenitorIDs)>0):
			#Get the location of where the MainBranchHostRootProgenitorIDs are
			MainBranchHostIndicator = np.in1d(AllRootProgenitors,MainBranchHostRootProgenitorIDs)
			#Set them so they point to the main branch
			AllBranchIndicators[MainBranchHostIndicator] = 0




	else:
		#Otherwise just set it equal to the TreeRootProgenitors
		AllRootProgenitors = TreeRootProgenitors
		AllBranchIndicators = ProgenBranchIndicator

	if(len(AllRootProgenitors)==1):
		raise Exception("This tree only has 1 branch so cannot be plotted")


	print("This tree has ",len(AllRootProgenitors),"branches")
	print("Done finding all RootProgenitors in",time.time() - start)

	return AllRootProgenitors, AllBranchIndicators


def MoveBranches(perIndx,branchIndicator,depthIndicator,Indexes,Indx,depth):

	#First find where the branches point to the perIndx
	subIndexes = np.sort(np.where(branchIndicator==perIndx)[0])[::-1]

	#Find where Indexes equals subIndexes
	delIndexes = np.where(np.in1d(Indexes,subIndexes))[0]

	#Delete where the index equals the sIndx in the array
	Indexes=np.delete(Indexes,delIndexes)

	#Find where the perIndex is in the Indexes array
	ploc = np.where(Indexes==perIndx)[0]

	#Now insert all the subIndexes next to the perIndex
	Indexes = np.insert(Indexes,ploc+1,subIndexes)

	#Lets find how deep the branches are
	Indicator = perIndx

	#Start at the depth of subbranches
	depth = 1

	#Loop until we reach a depth of 0
	while(Indicator>0):
		Indicator = branchIndicator[int(Indicator)]
		depth+=1
		if(depth>100):
			raise SystemExit("Reached a depth > 100, please check if the tree has been constructed correctly")

	#Set all these branches depth they were found at
	depthIndicator[subIndexes]=depth

	return Indx,Indexes,depthIndicator





def createPlotArrays(opt,plotOpt,treedata,SelIndex,outdir='',outputArrays=False):
	"""
	Function to generate the plot arrays

	"""

	plotData = {}

	AllRootProgenitors, AllBranchIndicators = getRootProgenitors(opt,treedata,SelIndex)

	if(len(AllRootProgenitors)==0):
		return [],[],[],[],[],[]


	start = time.time()
	print("Now creating the plot arrays")

	numBranches = np.int32(AllRootProgenitors.shape[0])

	#Intilize the plotting arrays for the data
	plotData["xposData"] = np.zeros([opt.Nsnaps,numBranches],dtype=np.float32) #The x-position of the point
	plotData["SizeData"] = np.zeros([opt.Nsnaps,numBranches],dtype=treedata["Snap_%03d" %opt.startSnap]["SizeData"].dtype) #The size of the data-point
	plotData["ColData"] = np.zeros([opt.Nsnaps,numBranches],dtype=treedata["Snap_%03d" %opt.startSnap]["ColData"].dtype) #The color of the point
	if(plotOpt.overplotdata):
		plotData["OverPlotData"] = np.zeros([opt.Nsnaps,numBranches],dtype=treedata["Snap_%03d" %opt.startSnap]["OverPlotData"].dtype)



	pos = np.zeros(3,dtype=np.float32)
	startpos = np.zeros(3,dtype=np.float32)
	prevpos = np.zeros(3,dtype=np.float32)

	# Get the snapshot which the main branch comes into existance
	MainBranchisnap = int(AllRootProgenitors[0]/opt.HALOIDVAL)

	#Store the main branches position and Radius
	mainBranchPos = np.zeros([opt.Nsnaps,3],dtype=np.float32)
	mainBranchRadius = np.zeros(opt.Nsnaps,dtype=np.float32)
	mainBranchIDs = np.zeros(opt.Nsnaps,dtype=np.int64)

	# print("ibranch snapKey haloID descProgenitor descendant branchRootTail descRootDescendant currRootDescendant descRootProgenitor currRootProgenitor")


	#Now we want to walk up from the AllRootProgenitors to build the tree and create the arrays
	for ibranch in range(numBranches):
		#Extract the Root Progenitor for this branch
		haloID = AllRootProgenitors[ibranch]
		
		#Move up the branch setting the plotting data as we go
		while(True):
			#Extract the information for this halo
			index = int(haloID%opt.HALOIDVAL-1)
			snap = int(haloID/opt.HALOIDVAL)


			snapKey = "Snap_%03d" %snap
			isnap = snap - opt.startSnap

			#Extract the position for this halo
			pos[:] =  treedata[snapKey]["Pos"][index]

			#If at the root progenitor set the previous position to the current
			if(haloID==AllRootProgenitors[ibranch]):
				prevpos[:]=pos

			for k in range(3):
				xposdiff = pos[k]-prevpos[k]
				if(abs(xposdiff)>0.8*opt.boxsize):
					pos[k]=pos[k]-opt.boxsize if(xposdiff>0) else pos[k]+opt.boxsize
			#Store its descendant
			descendant = treedata[snapKey]["Descendant"][index]

			#Set the size and col data to what was requested
			plotData["SizeData"][isnap,ibranch] = treedata[snapKey]["SizeData"][index]

			if(plotOpt.WWflag):
				plotData["ColData"][isnap,ibranch] = -2 if(treedata[snapKey]["WWHaloFlag"][index]) else treedata[snapKey]["ColData"][index]
			else:
				plotData["ColData"][isnap,ibranch] = treedata[snapKey]["ColData"][index]

			if(plotOpt.overplotdata):
					plotData["OverPlotData"][isnap,ibranch] = treedata[snapKey]["OverPlotData"][index]

			#If on the main branch
			if(ibranch==0):
				#Store both its position and Radius
				mainBranchPos[isnap,:] = pos
				mainBranchRadius[isnap] = cosFuncs.coshaloMvirToRvir(treedata[snapKey]["Mass"][index],Munit = opt.Munit,z= treedata[snapKey]["Redshift"]) 
				#If at the root progenitor then store its start position and set its position to be small (non-zero)
				if(haloID==AllRootProgenitors[ibranch]):
					startpos[:] = pos
					plotData["xposData"][isnap,ibranch] = 1e-5
				#Otherwise track its position relative to its start position
				else:
					plotData["xposData"][snap,ibranch] = np.sqrt((np.sum((pos - startpos))**2))

				#Store the main branches IDs
				mainBranchIDs[isnap] = haloID


				#If on the main branch check if we have reached the root head
				if(haloID==descendant): break


			else:

				#If the halo exist before the main branch set it to the limit of the plot
				if(snap<MainBranchisnap):
					plotData["xposData"][isnap,ibranch] = -1.0
				#Otherwise plot its distance relative to the main branch Radius at that snapshot
				else:

					if(mainBranchRadius[isnap]==0):
						print("The main branch Radius cannot be found for snapshot",snap,"it will be interpolated from the data for this snapshot")
						
						sel = mainBranchRadius>0
						snaps = np.where(sel)[0]
						f_Radius = interp1d(snaps,mainBranchRadius[sel])
						mainBranchRadius[isnap] = f_Radius(isnap)

	
					plotData["xposData"][isnap,ibranch] = np.sqrt(np.sum(((pos - mainBranchPos[isnap]))**2)) / mainBranchRadius[isnap]
					

				# Check if the descendant is equal to the haloID so we have reached the end of this branch
				if(haloID==descendant):
					# print("haloID==descendant",treedata[snapKey]["HaloID"][index],AllRootProgenitors[ibranch],treedata[snapKey]["RootDescendant"][index])
					
					#If we are not at the final snapshot then lets mark it as a branch which dies and doesn't merge with anything
					if(snap<(opt.Nsnaps-1)):
						AllBranchIndicators[ibranch]=-3



					break
				
				# Otherwise lets extract the descedants progenitor 
				descIndex = int(descendant%opt.HALOIDVAL-1)
				descSnap = int(descendant/opt.HALOIDVAL)
				descSnapKey = "Snap_%03d" %descSnap
				descProgenitor = treedata[descSnapKey]["Progenitor"][descIndex]

				#Now lets see if the progenitor points back to the haloID, if not then it has merged with another branch
				if(haloID!=descProgenitor): 

					#Lets find the RootProgenitor for this branch:
					branchRootTail = treedata[descSnapKey]["RootProgenitor"][descIndex]
					#Lets see if it doesn't point to the mainBranches RootProgenitor
					if(branchRootTail!=AllRootProgenitors[0]):
						#Otherwise lets track which branch in the RootProgenitor list it points to
						sel = np.where(AllRootProgenitors==branchRootTail)[0]


						
						if(sel==ibranch):
							#raise SystemExit("Something has gone wrong in the construction of the tree")
							break


						# If it merged with something that has aleady been found
						if(len(sel)==1):
							#Insert that into the branch indicators which one it points to
							AllBranchIndicators[ibranch]=sel
						else:
							AllBranchIndicators[ibranch]=-4

					break

			#Now move onto the next halo up the branch
			prevpos[:] = pos
			haloID = descendant

	##### Now lets sort the arrays using the data in the plotData["SizeData"] array

	#Store the maximum size in each branch
	BranchMaxSize = np.amax(plotData["SizeData"],axis=0)

	#Sort the branches so the largest branch are next to the main branch
	sortIndx = np.argsort(BranchMaxSize)[::-1].astype(np.int32)

	if(sortIndx[0]!=0):
		sel = sortIndx==0
		sortIndx[0],sortIndx[sel] = sortIndx[sel],sortIndx[0]


	#Now apply that sorting to the data
	depthIndicator = np.ones(numBranches,dtype=np.int32)
	depthIndicator[0]=0

	#Keep track  of which index we are at
	Indx = 0

	#Keep looping until we reach the end index
	while(Indx<(len(sortIndx)-1)):

		#Sort the branches so the progenitor branches are next to their merging branches
		Indx,sortIndx,depthIndicator = MoveBranches(sortIndx[Indx],AllBranchIndicators,depthIndicator,sortIndx,Indx+1,2)


	#Apply this sorting to the data
	plotData["xposData"] = plotData["xposData"][:,sortIndx]
	plotData["SizeData"] = plotData["SizeData"][:,sortIndx]
	plotData["ColData"] = plotData["ColData"][:,sortIndx]
	if(plotOpt.overplotdata):
		plotData["OverPlotData"] = plotData["OverPlotData"][:,sortIndx]
	AllBranchIndicators = AllBranchIndicators[sortIndx]
	depthIndicator = depthIndicator[sortIndx]



	#Set all the branches marked as subhalos to a depth of -1
	depthIndicator[AllBranchIndicators<-1]=-1

	print("Done creating the plot arrays in",time.time()-start)

	# if(outputArrays):

	# 	print("Outputing tree indx %i to disk" %SelIndex)
	# 	if(plotOpt.overplotdata):
	# 		np.savez_compressed(outdir+"/%i" %SelIndex, sortIndx=sortIndx, branchIndicators=AllBranchIndicators,depthIndicator=depthIndicator,mainBranchIDs=mainBranchIDs,xposData=plotData["xposData"],sizeData=plotData["SizeData"],colData=plotData["ColData"],OverPlotData=plotData["OverPlotData"])
	# 	else:
	# 		np.savez_compressed(outdir+"/%i" %SelIndex, sortIndx=sortIndx, branchIndicators=AllBranchIndicators,depthIndicator=depthIndicator,mainBranchIDs=mainBranchIDs,xposData=plotData["xposData"],sizeData=plotData["SizeData"],colData=plotData["ColData"],OverPlotData=plotData["OverPlotData"])

	# else:
	
	return plotData,AllBranchIndicators,depthIndicator,sortIndx,mainBranchIDs,mainBranchRadius



            
























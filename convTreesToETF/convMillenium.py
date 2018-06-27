import h5py
import numpy as np 
import time


def setProgenRootHeads(mainProgenSnap,mainHaloID,mainProgenIndex,BranchRootDesc,treedata,HALOIDVAL):

	mainProgenSnapKey = "Snap_%03d" %mainProgenSnap

	#Find which halos point to this halo as a descendant
	progenIndexes = np.where(halodata[mainProgenSnapKey]["Descendant"]==mainHaloID)[0]



	#Get the main branches other progenitors
	progenIndexes = progenIndexes[progenIndexes!=mainProgenIndex]


	#Loop over the halos progenitors
	for progenIndex in progenIndexes:


		# Set the progen snapshot
		iprogenSnap = progenSnap
		iprogenSnapKey = "Snap_%03d" %iprogenSnap

		#Set its EndDescendant and this branches depth
		treedata[iprogenSnapKey]["EndDescendant"][progenIndex] = BranchRootDesc


		#Loop until the get to the last progenitor
		while(haloID!=progenID):	

			#Get the progenitor's haloID and its progenitor
			progenID = halodata[iprogenSnapKey]["Progenitor"][int(progenIndex)]
			haloID = halodata[iprogenSnapKey]["HaloID"][int(progenIndex)]


			#Find the progenitors snapshot
			iprogenSnap = int(progenID/HALOIDVAL)
			iprogenSnapKey = "Snap_%03d" %iprogenSnap

			#Find the index where this progenitor's haloID is
			progenIndex = np.where(treedata[iprogenSnapKey]["HaloID"] == progenID)[0]


			#Set its EndDescendant and its depth
			treedata[iprogenSnapKey]["EndDescendant"][progenIndex] = BranchRootDesc

			#Check if halo has any progeintors and walk down its tree
			if(np.sum(halodata[iprogenSnapKey]["Descendant"]==haloID)>1):
				treedata = getProgenIndexes(opt,ibranch+branchOffset,iprogenSnap,haloID,progenIndex,treedata)
			

	return treedata


def convMilleniumToMTF(filename,numsnaps,snapKey,scalefactorKey,fieldsDict,HALOIDVAL = 10000000000000):


	halodata = loadMilleniumData(filename,numsnaps,snapKey,scalefactorKey,fieldsDict)

	start = time.time()

	numhalos = np.zeros(numsnaps,dtype=np.int64)

	requiredFields = ["HaloID","StartProgenitor","Progenitor","Descendant","EndDescendant","M200crit","Pos","HostHaloID","Redshift"]
	extraFields = [field for field in fieldsDict.keys() if field not in requiredFields] 

	treedata = {"Snap_%03d" %snap:{} for snap in range(numsnaps)}

	for snap in range(numsnaps):

		snapKey = "Snap_%03d" %snap

		numhalos[snap] = len(halodata[snapKey]["HaloID"])

		treedata[snapKey]["EndDescendant"] = np.zeros(numhalos[snap],dtype=np.int64)
		treedata[snapKey]["Descendant"] = np.zeros(numhalos[snap],dtype=np.int64)
		treedata[snapKey]["HaloID"] = np.zeros(numhalos[snap],dtype=np.int64)
		treedata[snapKey]["Progenitor"] = np.zeros(numhalos[snap],dtype=np.int64)
		treedata[snapKey]["StartProgenitor"] = np.zeros(numhalos[snap],dtype=np.int64)

		print(list(halodata[snapKey].keys()))
		treedata[snapKey]["Redshift"] = halodata[snapKey]["Redshift"]
		halodata[snapKey]["Redshift"] = None
		treedata[snapKey]["M200crit"] = halodata[snapKey]["M200crit"]
		halodata[snapKey]["M200crit"] = None
		treedata[snapKey]["Pos"] = halodata[snapKey]["Pos"]
		halodata[snapKey]["Pos"] = None
		for extraField in extraFields:
			treedata[snapKey][extraField]  = halodata[snapKey][extraField]
			halodata[snapKey][extraField] = None



	print("Setting Halos Progenitors, descendants and StartProgenitors")


	for snap in range(numsnaps):

		snapKey = "Snap_%03d" %snap
		print("Doing snap",snap)

		


		for ihalo in range(numhalos[snap]):

			if(treedata[snapKey]["HaloID"][ihalo]==0):

				#Set the new haloID
				treedata[snapKey]["HaloID"][ihalo] = snap*HALOIDVAL + ihalo + 1

				#Set the data for the current halo 
				currIndex = ihalo
				currSnap = snap
				currSnapKey = "Snap_%03d" %currSnap


				#Set the Branches RootTail
				BranchRootTail = treedata[currSnapKey]["HaloID"][currIndex]
				treedata[currSnapKey]["StartProgenitor"][currIndex] = BranchRootTail
				treedata[currSnapKey]["Progenitor"][currIndex]

				#Loop over this branch setting the Descendant, Progenitor and StartProgenitor
				while(True):

					#Extract the descendant
					DescID = halodata[currSnapKey]["Descendant"][currIndex]

					#Check to see if we have reached the EndDescendant of this branch
					if(DescID==-1):
						treedata[currSnapKey]["Descendant"][currIndex] = treedata[currSnapKey]["HaloID"][currIndex]
						treedata[currSnapKey]["EndDescendant"][currIndex] = treedata[currSnapKey]["HaloID"][currIndex]
						break

						
					#Get the descendants snapshot
					descSnap = int(DescID/HALOIDVAL)
					descSnapKey = "Snap_%03d" %descSnap


					#Get the index for the descedant in the next snapshot
					descIndex = np.where(halodata[descSnapKey]["HaloID"]==DescID)[0]

					#Check to see if this Progenitor ID has already be done
					if(treedata[descSnapKey]["Progenitor"][descIndex]>0):
						break

					#Convert the descendant ID to updated ID
					treedata[currSnapKey]["Descendant"][currIndex] = descSnap*HALOIDVAL + descIndex + 1

					#Set the halo ID for the descendant halo
					treedata[descSnapKey]["HaloID"][descIndex] = treedata[currSnapKey]["Descendant"][currIndex]

					#Set the tail and RootTail only if main progenitor or it involves a is a filler halo (has ID > 1e16)
					if(halodata[descSnapKey]["Progenitor"][descIndex]==halodata[currSnapKey]["HaloID"][currIndex]):
						treedata[descSnapKey]["StartProgenitor"][descIndex] = BranchRootTail
						treedata[descSnapKey]["Progenitor"][descIndex] = treedata[currSnapKey]["HaloID"][currIndex]

					#Move up so the current halo is the descendant
					currIndex = descIndex
					currSnap = descSnap
					currSnapKey = "Snap_%03d" %currSnap

	print("Done seeting the Progenitors,Descendants and StartProgenitors in",time.time()-start)

	print("Setting all the RootHead ID's")

	#Now we have built the Progenitors, Descendants and StartProgenitors, we need to now go back down the branches setting the RootDescedant for the halos
	for snap in range(numsnaps-1,-1,-1):

		print("Doing snap",snap)

		snapKey = "Snap_%03d" %snap


		for ihalo in range(numhalos[snap]):



			if(treedata[snapKey]["EndDescendant"][ihalo]):

				BranchRootDesc = treedata[snapKey]["EndDescendant"][ihalo] 

				mainProgenSnap = snap
				mainProgenSnapKey = "Snap_%03d" %mainProgenSnap
				mainProgenIndex = ihalo

				mainHaloID = treedata[mainProgenSnapKey]["HaloID"][mainProgenIndex]
				mainProgenID = treedata[mainProgenSnapKey]["Progenitor"][mainProgenIndex]

				while(mainHaloID!=mainProgenID):


					mainHaloID = treedata[mainProgenSnapKey]["HaloID"][mainProgenIndex]

					mainProgenID = treedata[mainProgenSnapKey]["Progenitor"][mainProgenIndex]

					mainProgenSnap = int(mainProgenID/HALOIDVAL)

					mainProgenSnapKey = "Snap_%03d" % mainProgenSnap

					mainProgenIndex = np.where(treedata[mainProgenSnapKey]["HaloID"]==mainProgenID)[0]

					treedata = setProgenRootHeads(mainProgenSnap,mainhaloID,mainProgenIndex,BranchRootDesc,treedata,HALOIDVAL)

					treedata[mainProgenSnapKey]["EndDescendant"][mainProgenIndex] = BranchRootDesc


	print("Done setting the EndDescendants in",time.time()-start)

	return treedata

def loadMilleniumData(filename,numsnaps,snapKey,scalefactorKey,dataKeys):

	hdffile = h5py.File(filename,"r")

	halodata = {"Snap_%03d" %snap:{} for snap in range(numsnaps)}

	# Load in the snapshot date from the catalogue
	snapData = hdffile[snapKey][:]

	scalefactorData = hdffile[scalefactorKey][:]

	pad = np.zeros(numsnaps - len(scalefactorData))

	scalefactorData = np.concatenate([pad,scalefactorData])

	#Load in all the requested datasets for fast access
	tmpData = {key:hdffile[key][:]   for key in dataKeys.values()}


	hdffile.close()

	# Now we can load in the data, splitting it across snapshots for quicker searching
	for snap in range(numsnaps):
		print("Loading snap",snap)

		snapKey = "Snap_%03d" %snap

		snapSel = np.where(snapData==snap)[0].tolist()


		halodata[snapKey]["Redshift"] = scalefactorData[snap]

		for dataKey in dataKeys.keys():
			halodata[snapKey][dataKey] = tmpData[dataKeys[dataKey]][snapSel]




	
	

	
	return halodata



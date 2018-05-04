import numpy as np
import sys
import h5py
import time
import multiprocessing as mp
import array
from functools import partial


def convAHFToMTF(startSnap,endSnap,haloFilelist,treeFilelist,fieldsDict):

	numsnaps = endSnap - startSnap +1

	Redshift,halodata = ReadInHaloFilesAcrossSnapshots(startSnap,endSnap,haloFilelist,fieldsDict)

	treedata = ReadHaloMergerTreeAcrossSnapshots(startSnap,endSnap,treeFilelist)

	halodata = convToMTF(startSnap,endSnap,halodata,treedata)

	return Redshift,halodata


def ReadInHaloFilesAcrossSnapshots(startSnap,endSnap,filelist,fieldsDict):
	"""
	Function to read in the AHF halo data across many files

	"""

	numsnaps = abs(endSnap - startSnap + 1)

	fields = list(fieldsDict.values())

	snapfilelist=open(filelist,"r")
	start=time.clock()  
	halodata={}


	filename=snapfilelist.readline().strip()
	halofile = open(filename,"r")
	props = halofile.readline().strip("#").split()

	halofile.close()

	allfieldnames = []
	allfielddtypes = []
	for field in fields:

		if(field[0]!=""):
			if("," not in field[0]):
				allfieldnames.append(field[0])
				allfielddtypes.append((field[0],field[1]))

			else:

				for ifield in field[0].split(","):
					allfieldnames.append(ifield)
					allfielddtypes.append((ifield,field[1]))


	colReadIndx =  np.where(np.in1d(props,allfieldnames))[0]

	fieldnames = [props[indx] for indx in colReadIndx]
	fielddtypes = [allfielddtypes[allfieldnames.index(field)]  for field in fieldnames]

	snapfilelist.close()

	snapfilelist=open(filelist,"r")

	Redshift = np.zeros(endSnap-startSnap+1,dtype=int)

	for snap in range(endSnap,startSnap-1,-1):

		snapKey = "Snap_%03d" %snap

		halodata[snapKey] = {}

		filename=snapfilelist.readline().strip()
		print("Reading ",filename)

		isnap =  snap - startSnap
		# Extract the radshift from the filename
		Redshift[isnap] = float(filename[filename.find("z")+1:filename.find(".AHF")])

		halofile = open(filename,"r")

		with np.warnings.catch_warnings():
			np.warnings.filterwarnings('ignore')
			tmpData = np.loadtxt(halofile,usecols = colReadIndx,dtype = fielddtypes,ndmin=2 )


		for MTFfieldname in fieldsDict.keys():

			fieldname = fieldsDict[MTFfieldname][0]

			if(fieldname!=""):
				if("," not in fieldname):
					halodata[snapKey][MTFfieldname] = tmpData[fieldname].flatten()

				elif(MTFfieldname=="Pos"):
					splitfields  = fieldname.split(",")
					halodata[snapKey][MTFfieldname]=np.column_stack([tmpData[splitfield] for splitfield in splitfields])/1000

				else:
					splitfields  = fieldname.split(",")
					halodata[snapKey][MTFfieldname]=np.column_stack([tmpData[splitfield] for splitfield in splitfields])
		halofile.close()

	print("Halo data read in ",time.clock()-start)

	snapfilelist.close()

	return Redshift,halodata


def ReadHaloMergerTreeAcrossSnapshots(startSnap,endSnap,filelist):
	"""
	Function to read the Merger Tree from AHF across snapshots
	"""
	start=time.clock()  
	snapfilelist=open(filelist,"r")
	# halodata = 

	treedata = {"Snap_%03d" %startSnap:{"HaloID": [], "Num_progen": [], "Progenitors": [],"SharedNpart":[],"mainProgenitor":[]}}
	for snap in range(endSnap,startSnap,-1):
		filename=snapfilelist.readline().strip()
		print("Reading treefile",filename)
		treefile=open(filename,"r")
		header1=treefile.readline()
		header2=treefile.readline()
		treefile_idx = open(filename+"_idx","r")
		header_idx = treefile_idx.readline()
	   
		snapKey = "Snap_%03d" %snap

		treedata[snapKey] = {"HaloID": [], "Num_progen": [], "Progenitors": [],"SharedNpart":[],"mainProgenitor":[]}
		 
		j=0

		while(True):
			try:
				HostID,HostNpart,Num_progen=treefile.readline().strip().split("  ")
				HostID,Progenitor = treefile_idx.readline().strip().split(" ")
			except ValueError:
				break
			treedata[snapKey]["HaloID"].append(np.int64(HostID))
			treedata[snapKey]["Num_progen"].append(np.int(Num_progen))
			treedata[snapKey]["mainProgenitor"].append(np.int64(Progenitor))
			treedata[snapKey]["Progenitors"].append([])
			treedata[snapKey]["SharedNpart"].append([])
			for iprogen in range(int(Num_progen)):
				line=treefile.readline().strip().split("  ")
				[SharedNpart,SatID,SatNpart]=line
				treedata[snapKey]["Progenitors"][j].append(np.int64(SatID))
				treedata[snapKey]["SharedNpart"][j].append(np.int64(SharedNpart))
			treedata[snapKey]["Progenitors"][j] = np.array(treedata[snapKey]["Progenitors"][j])
			treedata[snapKey]["SharedNpart"][j] = np.array(treedata[snapKey]["SharedNpart"][j])
			j+=1
		treedata[snapKey]["HaloID"]=np.array(treedata[snapKey]["HaloID"])
		treedata[snapKey]["Num_progen"]=np.array(treedata[snapKey]["Num_progen"])

	print("Tree data read in ",time.clock()-start)
	return treedata	


def walkDownProgenBranches(snap,halodata,treedata,MTFdata,Progenitors,Descendant,TreeRootDescendant,HALOIDVAL,startSnap,depth,treeProgenIndex=0):

	#print(Descendant,Progenitors,treeProgenIndex)

	

	

	for haloID in Progenitors:

		haloSnap = snap
		haloSnapKey = "Snap_%03d" %snap 

		# Lets set these progen branches to point up to the descendant and RootDescendant halo
		haloIndex = np.where(halodata[haloSnapKey]["HaloID"]==haloID)[0].astype(int)
		
		if(MTFdata[haloSnapKey]["HaloID"][haloIndex]==0):

			MTFdata[haloSnapKey]["Descendant"][haloIndex] = Descendant

			MTFdata[haloSnapKey]["RootDescendant"][haloIndex] = TreeRootDescendant

			while(True):

				#Lets set the ID for this halo
				MTFdata[haloSnapKey]["HaloID"][haloIndex] = haloSnap * HALOIDVAL + haloIndex +1
				#if(depth==0): print("Doing halo",MTFdata[haloSnapKey]["HaloID"][haloIndex],haloID)

				#If we are at the endsnap of the simulation then we don't have to search for the final progenitor
				if(haloSnap==startSnap):
					MTFdata[haloSnapKey]["Progenitor"][haloIndex] = MTFdata[haloSnapKey]["HaloID"][haloIndex]
					break

				# Find where it exists in the treedata
				treeProgenIndx = np.where(treedata[haloSnapKey]["HaloID"]==haloID)[0]

				# If it doesn't exist in the treedata lets break out of the loop
				if(treeProgenIndx.size==0):
					#print("Reached RootProgenitor for halo",MTFdata[haloSnapKey]["HaloID"][haloIndex])
					MTFdata[haloSnapKey]["Progenitor"][haloIndex] = MTFdata[haloSnapKey]["HaloID"][haloIndex]
					break

				# Lets find all the progenitors for this halo
				Progenitors = treedata[haloSnapKey]["Progenitors"][int(treeProgenIndx)]

				#Lets find its main progenitor
				mainProgenitor = treedata[haloSnapKey]["mainProgenitor"][int(treeProgenIndx)]

				# Make sure the main progenitor isn't in the list of progenitors
				SecondaryProgenitors = Progenitors[Progenitors!=mainProgenitor]

				MTFdata = walkDownProgenBranches(haloSnap - 1,halodata,treedata,MTFdata,SecondaryProgenitors,MTFdata[haloSnapKey]["HaloID"][haloIndex],TreeRootDescendant,HALOIDVAL,startSnap,depth+1,treeProgenIndx)

				# Lets find the info for its progenitor
				progenSnap = haloSnap - 1
				progenSnapKey = "Snap_%03d" %progenSnap
				progenIndex = int(np.where(halodata[progenSnapKey]["HaloID"]==mainProgenitor)[0])

				# Ltes set the current halos progenitor as the main progenitor
				MTFdata[haloSnapKey]["Progenitor"][haloIndex] = progenSnap * HALOIDVAL + progenIndex +1

				# Lets set the current 
				MTFdata[progenSnapKey]["Descendant"][progenIndex] = MTFdata[haloSnapKey]["HaloID"][haloIndex]
				MTFdata[progenSnapKey]["RootDescendant"][progenIndex] = TreeRootDescendant

				#Now lets move to the main progenitor halo
				haloSnap = progenSnap
				haloSnapKey = progenSnapKey
				haloIndex = progenIndex
				haloID= mainProgenitor

	return MTFdata

def SetProgenitorandDescendants(snap,startSnap,numhalos,halodata,treedata,MTFdata,HALOIDVAL):

	progenSnap = snap - 1

	isnap =  snap - startSnap
	print("Doing snap",snap,isnap)



	snapKey = "Snap_%03d" %snap
	progenSnapKey = "Snap_%03d" %(snap-1)


	snapData = MTFdata[snapKey]
	progenSnapData = MTFdata[progenSnapKey]

	

	for ihalo in range(numhalos[isnap]):
		

		HaloID = halodata[snapKey]["HaloID"][ihalo]

		HostID = halodata[snapKey]["HostHaloID"][ihalo]
		if(HostID!=0):

			#Find the location of the host in the catalogue
			HostLoc = int(np.where(halodata[snapKey]["HaloID"]==HostID)[0])

			#Adjust the hosts ID
			snapData["HostHaloID"][ihalo] = snap * HALOIDVAL + HostLoc + 1

		#Find where in the treedata the haloID equal this haloID
		treeProgenIndx =np.where(treedata[snapKey]["HaloID"]==HaloID)[0]

		#If it does not exist in the treedata then that branch no longer exists
		if(treeProgenIndx.size==0):
			#Set its progenitor to point back to itself
			snapData["Progenitor"][ihalo] = snapData["HaloID"][ihalo]
			continue
		

		# Find all this halo progenitors
		Progenitors = treedata[snapKey]["Progenitors"][int(treeProgenIndx)]

		# Find the main progenitor for this halo
		mainProgenitor = treedata[snapKey]["mainProgenitor"][int(treeProgenIndx)]
		# print("This haloID",haloID, " mainProgenitor is", mainProgenitor)

		# Remove the main progenitor from the list of progenitors
		Progenitors = Progenitors[Progenitors!=mainProgenitor]

		# Set the data for the main progenitor
		progenIndex = int(np.where(halodata[progenSnapKey]["HaloID"]==mainProgenitor)[0])
		snapData["Progenitor"][ihalo] = progenSnap * HALOIDVAL + progenIndex +1

		#Lets set the descendant of the progenitor halo
		progenSnapData["Descendant"][progenIndex] = snapData["HaloID"][ihalo] 

		#Now lets set the rest of the progenitors descendants to be this halo
		for progen in Progenitors:

			progenIndex = int(np.where(halodata[progenSnapKey]["HaloID"]==progen)[0])

			progenSnapData["Descendant"][progenIndex] = snapData["HaloID"][ihalo] 

	
	MTFdata[snapKey] = snapData


	deadSel = progenSnapData["Descendant"]==0
	progenSnapData["Descendant"][deadSel]=progenSnapData["HaloID"][deadSel]
	MTFdata[progenSnapKey] = progenSnapData
		

		
def SetProgenandDesc(ihalo,snap,startSnap,halodata,treedata,MTFdata,HALOIDVAL):

	
	snapKey = "Snap_%03d" %snap

	HostID = halodata[snapKey]["HostHaloID"][ihalo]
	if(HostID!=0):

		#Find the location of the host in the catalogue
		HostLoc = int(np.where(halodata[snapKey]["HaloID"]==HostID)[0])

		#Adjust the hosts ID
		MTFdata[snapKey]["HostHaloID"][ihalo] = snap * HALOIDVAL + HostLoc + 1

	if(MTFdata[snapKey]["HaloID"][ihalo]==0):
		


		MTFdata[snapKey]["HaloID"][ihalo] = snap * HALOIDVAL + ihalo +1
	

		TreeRootDescendant = MTFdata[snapKey]["HaloID"][ihalo]
		Descendant = TreeRootDescendant

		# print("Doing tree with RootDescendant", TreeRootDescendant,"of",numhalos[isnap])
		
		#Set the data for the current halo 
		haloIndex = ihalo
		haloSnap = snap
		haloSnapKey = snapKey
		haloID = halodata[haloSnapKey]["HaloID"][ihalo]


		while(True):

			#Set the haloID, descendant and the RootDescendant for this halo
			MTFdata[haloSnapKey]["Descendant"][haloIndex] = Descendant
			MTFdata[haloSnapKey]["RootDescendant"][haloIndex] = TreeRootDescendant

			# print("On haloID",MTFdata[haloSnapKey]["HaloID"][haloIndex])


			#If we are at the endsnap of the simulation then we don't have to search for the final progenitor
			if(haloSnap==startSnap):
				MTFdata[haloSnapKey]["Progenitor"][haloIndex] = MTFdata[haloSnapKey]["HaloID"][haloIndex]
				break


			#Find where in the treedata the haloID equal this haloID
			treeProgenIndx = np.where(treedata[haloSnapKey]["HaloID"]==haloID)[0]

			#If it does not exist in the treedata then that branch no longer exists
			if(treeProgenIndx.size==0):
				#Set its progenitor to point back to itself
				MTFdata[haloSnapKey]["Progenitor"][haloIndex] = MTFdata[haloSnapKey]["HaloID"][haloIndex]
				break


			# Find all this halo progenitors
			Progenitors = treedata[haloSnapKey]["Progenitors"][int(treeProgenIndx)]

			# Find the main progenitor for this halo
			mainProgenitor = treedata[haloSnapKey]["mainProgenitor"][int(treeProgenIndx)]
			# print("This haloID",haloID, " mainProgenitor is", mainProgenitor)

			# Remove the main progenitor from the list of progenitors
			SecondaryProgenitors = Progenitors[Progenitors!=mainProgenitor]

			# Walk down the progenitor branches setting the Progenitor, Descedant and RootDesendant for each halo
			MTFdata = walkDownProgenBranches(haloSnap - 1,halodata,treedata,MTFdata,SecondaryProgenitors,MTFdata[haloSnapKey]["HaloID"][haloIndex],TreeRootDescendant,HALOIDVAL,startSnap,0)

			# Set the data for the main progenitor
			progenSnap = haloSnap - 1
			progenSnapKey = "Snap_%03d" %progenSnap
			progenIndex = int(np.where(halodata[progenSnapKey]["HaloID"]==mainProgenitor)[0])

			MTFdata[haloSnapKey]["Progenitor"][haloIndex] = progenSnap * HALOIDVAL + progenIndex +1
			Descendant = MTFdata[haloSnapKey]["HaloID"][haloIndex]
	
			haloSnap = progenSnap
			haloSnapKey = progenSnapKey
			haloIndex = progenIndex
			haloID= mainProgenitor

			MTFdata[haloSnapKey]["HaloID"][haloIndex] = haloSnap * HALOIDVAL + haloIndex +1

		# print("Done tree with RootDescendant",TreeRootDescendant,"in",time.time()-starthalo)




def SetProgenandDescParallel(snap,halochunk,startSnap,endSnap,halodata,treedata,MTFdata,mpMTFdata,lock,HALOIDVAL):

	name = mp.current_process().name
	print(name,"is doing",np.min(halochunk),"to",np.max(halochunk))

	for ihalo in halochunk:
		SetProgenandDesc(ihalo,snap,startSnap,halodata,treedata,MTFdata,HALOIDVAL)

	print(name,"is on to copying the data")

	#Aquire the lock for the mpMTFdata
	lock.acquire()

	#Lets copy the local MTFdata to the global mpMTFdata
	for snap in range(startSnap,endSnap+1):
		snapKey = "Snap_%03d" %snap

		#Find the additions to the dataset
		selIndexes = np.where((MTFdata[snapKey]["HaloID"]>0))

		for key in MTFdata[snapKey].keys():

			#Find the typecode for the datatype
			dtype = np.sctype2char(MTFdata[snapKey][key])

			#Get the data currently stored in mpMTFdata
			tmpData = np.frombuffer(mpMTFdata[snapKey][key][:],dtype=dtype)

			#Put these changes into the tmp data
			tmpData[selIndexes] = MTFdata[snapKey][key][selIndexes]

			#Insert this array back into the mpMTFdata
			mpMTFdata[snapKey][key][:] = array.array(dtype,tmpData)

	#Release the lock on the data
	lock.release()

	#Delete this processes local copy of the MTFdata
	del MTFdata

	print(name,"is done")


def convToMTF(startSnap,endSnap,halodata,treedata,HALOIDVAL = 1000000000000):




	totstart = time.time()

	requiredFields = ["HaloID","RootProgenitor","Progenitor","Descendant","RootDescendant","M200crit","Pos","HostHaloID"]
	extraFields = [field for field in halodata["Snap_%03d" %startSnap].keys() if field not in requiredFields] 

	numsnaps = endSnap - startSnap + 1

	numhalos = np.zeros(numsnaps,dtype=np.int64)


	#Setup a dictionary for the MTF data
	MTFdata = {}


	for snap in range(endSnap,startSnap-1,-1):

		snapKey = "Snap_%03d" %snap
		isnap = endSnap - snap 

		# Setup another dictionary inside each snapshto to store the data
		MTFdata[snapKey] = {}

		#Find the number of halos at this snapshot
		numhalos[isnap] = halodata[snapKey]["HaloID"].size

		# Intialize arrays for the tree properties
		MTFdata[snapKey]["RootDescendant"] = np.zeros(numhalos[isnap],dtype=np.int64)
		MTFdata[snapKey]["Descendant"] = np.zeros(numhalos[isnap],dtype=np.int64)
		MTFdata[snapKey]["HaloID"] = np.zeros(numhalos[isnap],dtype=np.int64)
		MTFdata[snapKey]["Progenitor"] = np.zeros(numhalos[isnap],dtype=np.int64)
		MTFdata[snapKey]["RootProgenitor"] = np.zeros(numhalos[isnap],dtype=np.int64)
		MTFdata[snapKey]["HostHaloID"] = -1 * np.ones(numhalos[isnap],dtype=np.int64)


	print("Setting Progenitors, Descendant and RootDescendants for the branches")

	manager = mp.Manager()

	lock = manager.Lock()

	mpMTFdata=manager.dict({"Snap_%03d" %snap:manager.dict({key:manager.Array(np.sctype2char(MTFdata["Snap_%03d" %snap][key]),MTFdata["Snap_%03d" %snap][key]) for key in MTFdata["Snap_%03d" %snap].keys()}) for snap in range(startSnap,endSnap+1)})

	chunksize=10000

	for snap in range(endSnap,startSnap-1,-1):

		start = time.time()

		isnap = endSnap - snap

		snapKey = "Snap_%03d" %snap

		selUnset = MTFdata[snapKey]["HaloID"]==0

		ihalos = np.where(selUnset | (halodata[snapKey]["HostHaloID"]!=0))[0] 

		numhalounset = np.sum(selUnset)

		inumhalos = len(ihalos)

		print("Doing snap", snap,"with",inumhalos,"unset halos")


		if(numhalounset>2*chunksize):

			nthreads=int(min(mp.cpu_count(),np.ceil(inumhalos/float(chunksize))))
			nchunks=int(np.ceil(inumhalos/float(chunksize)/float(nthreads)))
			print("Using", nthreads,"threads to parse ",inumhalos," halos in ",nchunks,"chunks, each of size", chunksize)
			#now for each chunk run a set of proceses
			for j in range(nchunks):
				offset=j*nthreads*chunksize
				#if last chunk then must adjust nthreads
				if (j==nchunks-1):
					nthreads=int(np.ceil((inumhalos-offset)/float(chunksize)))

				
				#adjust last chunk
				if (j==nchunks-1):
					halochunk=[ihalos[offset+k*chunksize:offset+(k+1)*chunksize] for k in range(nthreads-1)]
					halochunk.append(ihalos[offset+(nthreads-1)*chunksize:inumhalos])
				else:
					halochunk=[ihalos[offset+k*chunksize:offset+(k+1)*chunksize] for k in range(nthreads)]
				#adjust last chunk
				if (j==nchunks-1):
					halochunk[-1]=range(offset+(nthreads-1)*chunksize,numhalos[isnap])
				#when calling a process pass not just a work queue but the pointers to where data should be stored
				processes=[mp.Process(target=SetProgenandDescParallel,args=(snap,halochunk[k],startSnap,endSnap,halodata,treedata,MTFdata,mpMTFdata,lock,HALOIDVAL)) for k in range(nthreads)]
				count=0
				for p in processes:
					p.start()
					count+=1
				for p in processes:
					#join thread and see if still active
					p.join()

			for snap in range(startSnap,endSnap+1):
				snapKey = "Snap_%03d" %snap
				# sel = np.where(MTFdata[snapKey]["HaloID"]>0)[0]
				for key in MTFdata[snapKey].keys():
					MTFdata[snapKey][key][:] = np.array(mpMTFdata[snapKey][key][:])


		else:

			for ihalo in ihalos:

				SetProgenandDesc(ihalo,snap,startSnap,halodata,treedata,MTFdata,HALOIDVAL)


		print("Done snap in",time.time()-start,np.sum(MTFdata["Snap_%03d" %snap]["HaloID"]==0))	

	print("Setting RootProgenitors")

	#Now we have set the Progenitors, Descendants and RootDescendants for all the halos we now need to
	#Walk back up the tree setting the RootProgenitors


	for snap in range(startSnap,endSnap+1):

		snapKey = "Snap_%03d" %snap
		isnap = endSnap - snap

		for ihalo in range(numhalos[isnap]):

			#First lets check if its root progenitor has been set
			if(MTFdata[snapKey]["RootProgenitor"][ihalo]==0):

				# If not extract the halo ID and set it as the root progenitor
				HaloID = MTFdata[snapKey]["HaloID"][ihalo]
				branchRootProgenitor = HaloID

				#Setting itself as it own root progenitor
				MTFdata[snapKey]["RootProgenitor"][ihalo] = branchRootProgenitor

				# Extract the halos desendant 
				Descendant = MTFdata[snapKey]["Descendant"][ihalo]
				descSnap = int(Descendant/HALOIDVAL)
				descSnapKey = "Snap_%03d" %descSnap
				descIndex = int(Descendant%HALOIDVAL-1)

				# Get the descendants progenitor
				DescendantProgen = MTFdata[descSnapKey]["Progenitor"][descIndex]

				#Lets check we haven't reached the end of the branch or if the halo is the main progenitor for the descendant halo
				while((HaloID!=Descendant) & (DescendantProgen==HaloID)):

					#Lets move to the descendant halo and set its root progenitor
					HaloID = Descendant
					haloSnap = descSnap
					haloSnapKey = descSnapKey
					haloIndex = descIndex
					MTFdata[haloSnapKey]["RootProgenitor"][haloIndex] = branchRootProgenitor

					#Extract the halos desendant and the desendants progenitor
					Descendant = MTFdata[haloSnapKey]["Descendant"][haloIndex]
					descSnap = int(Descendant/HALOIDVAL)
					descSnapKey = "Snap_%03d" %descSnap
					descIndex = int(Descendant%HALOIDVAL-1)
					DescendantProgen = MTFdata[descSnapKey]["Progenitor"][descIndex]					




	# Lets set the data in the extra fields
	for snap in range(endSnap,startSnap-1,-1):

		snapKey = "Snap_%03d" %snap

		MTFdata[snapKey]["Mass"] = halodata[snapKey]["Mass"]
		MTFdata[snapKey]["Rvir"] = halodata[snapKey]["Rvir"]
		MTFdata[snapKey]["Pos"] = halodata[snapKey]["Pos"]
		for extraField in extraFields:
				MTFdata[snapKey][extraField]  = halodata[snapKey][extraField]


	print("Done conversion in",time.time()-totstart)


	return MTFdata

				
















import numpy as np
import time
import multiprocessing as mp
import array


def convRockstarToMTF(snapColName,numProgName,redshiftColName,startSnap,endSnap,treefilelist,fieldsDict):

	Redshift,MTFdata = LoadRockstarIntoMTF(snapColName,numProgName,redshiftColName,startSnap,endSnap,treefilelist,fieldsDict)

	MTFdata = convToMTF(startSnap,endSnap,fieldsDict,MTFdata)


	return Redshift,MTFdata

def _getconv(dtype):
	""" Find the correct dtype converter. Adapted from matplotlib """



	typ = dtype.type
	if issubclass(typ, np.bool_):
		return lambda x: bool(int(x))
	if issubclass(typ, np.uint64):
		return np.uint64
	if issubclass(typ, np.int64):
		return np.int64
	if issubclass(typ, np.int32):
		return np.int32
	elif issubclass(typ, np.longdouble):
		return np.longdouble
	elif issubclass(typ, np.floating):
		return np.float32
	else:
		raise SystemExit("Incorrect data type")

def LoadRockstarIntoMTF(snapColName,numProgName,redshiftColName,startSnap,endSnap,treefilelist,fieldsDict,HALOIDVAL=1000000000000):

	treeFields = ["HaloID","RootProgenitor","Progenitor","Descendant","RootDescendant"]
	otherFields = [field for field in fieldsDict.keys() if field not in treeFields] 

	start = time.time()
	treelist = open(treefilelist,"r")

	filename=treelist.readline().strip()
	treefile = open(filename,"r")
	props = treefile.readline().strip("#").split()
	treefile.close()

	allfieldnames = [snapColName,numProgName,redshiftColName]
	allfielddtypes = [np.int32,np.int32,np.float32]
	for field in fieldsDict.values():

		if(field[0]!=""):
			if("," not in field[0]):
				allfieldnames.append(field[0])
				allfielddtypes.append(_getconv(np.dtype(field[1])))

			else:

				for ifield in field[0].split(","):
					allfieldnames.append(ifield)
					allfielddtypes.append(_getconv(np.dtype(field[1])))


	colReadIndx =  np.where(np.in1d(props,allfieldnames))[0]

	fieldnames = [props[indx] for indx in colReadIndx]
	fielddtypes = [allfielddtypes[allfieldnames.index(field)]  for field in fieldnames]

	treelist.close()

	treelist = open(treefilelist,"r")

	tmphalodata = {}

	MTFfieldnames = list(fieldsDict.keys()) + ["RootDescendant","RootProgenitor"]

	MTFdata = {"Snap_%03d" %snap:{field:[] for field in MTFfieldnames} for snap in range(startSnap,endSnap+1)}

	snapsIndex = np.zeros(endSnap-startSnap,dtype=int)

	Redshift = np.zeros(endSnap-startSnap+1 ,dtype=np.float32)

	for line in treelist:

		treefile = line.strip()

		print("Reading ",treefile)

		fp = open(treefile,"r")

		#Lets skip over the header and the number of trees in the file
		l = fp.readline()
		while(l[0]=="#"):
			l = fp.readline()

		numTrees=int(l)
		# Skip over the first tree number
		fp.readline()

		for i in range(numTrees):

			if(i%1000==0): print("On tree",i,"out of",numTrees)
			line = fp.readline()

			ihalo = 0
			prevNumProg=np.zeros([endSnap-startSnap],dtype=np.int64)

			if not line:
				break


			while(line[0]!="#"):
				line = line.strip().split()

				data = {fieldnames[j]:fielddtypes[j](line[indx]) for j,indx in enumerate(colReadIndx)}


				mainsnap = data[snapColName]
				mainsnapKey = "Snap_%03d" %mainsnap

				if(Redshift[mainsnap]==0):
					Redshift[mainsnap] = 1/data[redshiftColName] -1

				HaloID = mainsnap * HALOIDVAL + len(MTFdata["Snap_%03d" %mainsnap]["HaloID"]) + 1

				if(ihalo==0):
					TreeRootDescendant = HaloID
					Descendant = HaloID
					MTFdata[mainsnapKey]["Descendant"].append(Descendant)
					MTFdata[mainsnapKey]["RootDescendant"].append(TreeRootDescendant)

				MTFdata[mainsnapKey]["HaloID"].append(HaloID)
								
				
				for field in otherFields:
						 
					if("," not in fieldsDict[field][0]):
						MTFdata[mainsnapKey][field].append(data[fieldsDict[field][0]])
					else:
						MTFdata[mainsnapKey][field].append([data[ifield] for ifield in fieldsDict[field][0].split(",")])



				numProg = data[numProgName]

				if(numProg==0):
					MTFdata[mainsnapKey]["Progenitor"].append(HaloID)
					MTFdata[mainsnapKey]["RootProgenitor"].append(0)


				else:
					progenSnap = mainsnap - 1
					progenSnapKey = "Snap_%03d" %progenSnap

					MTFdata[mainsnapKey]["Progenitor"].append(progenSnap * HALOIDVAL + len(MTFdata[progenSnapKey]["HaloID"]) + prevNumProg[progenSnap]  + 1)
					MTFdata[mainsnapKey]["RootProgenitor"].append(0)

					MTFdata[progenSnapKey]["Descendant"].extend([HaloID]*numProg)
					MTFdata[progenSnapKey]["RootDescendant"].extend([TreeRootDescendant]*numProg)

				#move onto the next line
				line = fp.readline()

				if not line:
					break

				ihalo+=1
				prevNumProg[progenSnap]+=numProg

	# print(MTFdata["Snap_010"]["origID"])

	#Convert everything into array for easy indexing
	for snap in range(startSnap,endSnap+1):
		snapKey = "Snap_%03d" %snap
		for field in MTFdata[snapKey].keys():
			MTFdata[snapKey][field] = np.asarray(MTFdata[snapKey][field])

	# print(MTFdata["Snap_010"]["origID"])
	print("Done loading the data into EFT format in",time.time()-start)


	return Redshift,MTFdata


def SetHostID(ihaloset,snap,HostHaloID,origID,HALOIDVAL):
	
	for ihalo in ihaloset:

		#Lets Change the HostHaloID to the ETF halo ID
		if(HostHaloID[ihalo]>-1):

			#Find where the ID is in the original ID
			hostIndx = np.where(origID==HostHaloID[ihalo])[0]

			#Set the host ID
			HostHaloID[ihalo] = snap * HALOIDVAL + hostIndx + 1

	return HostHaloID


def SetHostIDParallel(ihaloset,snap,HostHaloID,origID,mpMTFdata,HALOIDVAL):


	name = mp.current_process().name
	print(name,"is doing",snap)

	HostHaloID=SetHostID(ihaloset,snap,HostHaloID,origID,HALOIDVAL)


	print(name,"is on to copying the data")

	snapKey = "Snap_%03d" %snap

	#Extract the data that is there
	tmpData = np.frombuffer(mpMTFdata[snapKey][:],dtype=np.int64)

	#Put these changes into the tmp data
	tmpData[selIndexes] = HostHaloID[ihaloset]

	#Insert this array back into the mpMTFdata
	mpMTFdata[snapKey][:] = array.array("l",tmpData)

	#Delete this processes local copy of the MTFdata
	del MTFdata

	print(name,"is done")




def convToMTF(startSnap,endSnap,fieldsDict,MTFdata,HALOIDVAL=1000000000000):


	totstart = time.time()

	
	print("Setting RootProgenitors")

	#Now we have set the Progenitors, Descendants and RootDescendants for all the halos we now need to
	#Walk back up the tree setting the RootProgenitors


	for snap in range(startSnap,endSnap+1):

		start = time.time()

		snapKey = "Snap_%03d" %snap

		numhalos = len(MTFdata[snapKey]["HaloID"])

		for ihalo in range(numhalos):

			HostID = MTFdata[snapKey]["HostHaloID"][ihalo]
			if(HostID>-1):

				#Find the location of the host in the catalogue
				HostLoc = int(np.where(MTFdata[snapKey]["origID"]==HostID)[0])

				#Adjust the hosts ID
				MTFdata[snapKey]["HostHaloID"][ihalo] = snap * HALOIDVAL + HostLoc + 1


			#First lets check if its root progenitor has been set
			if(MTFdata[snapKey]["RootProgenitor"][ihalo]==0):

				# If not extract the halo ID and set it as the root progenitor
				HaloID = MTFdata[snapKey]["HaloID"][ihalo]
				branchRootProgenitor = HaloID

				MTFdata[snapKey]["RootProgenitor"][ihalo]=branchRootProgenitor

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

		print("Done snap",snap,"in",time.time()-start)



	print("Done setting RootProgenitors in",time.time()-totstart)



	# print("Now setting the HostHaloID")


	# chunksize=500000

	# isnap = startSnap
	
	# ihaloset = np.where(MTFdata[snapKey]["HostHaloID"]>-1)[0] 

	# while(len(ihaloset)<chunksize):

	# 	print("doing snap",isnap)
	# 	MTFdata[snapKey]["HostHaloID"] = SetHostID(ihaloset,isnap,MTFdata[snapKey]["HostHaloID"],MTFdata[snapKey]["origID"],HALOIDVAL)

	# 	isnap+=1
	# 	ihaloset = np.where(MTFdata[snapKey]["HostHaloID"]>-1)[0] 


	# if(isnap<endSnap):

	# 	manager = mp.Manager()

	# 	mpMTFdata=manager.dict({"Snap_%03d" %snap:manager.Array(np.sctype2char(MTFdata["Snap_%03d" %snap]["HostHaloID"]),MTFdata["Snap_%03d" %snap]["HostHaloID"]) for snap in range(isnap,endSnap+1)})

		
	# 	inumsnaps = endSnap-isnap + 1

	# 	nthreads=min(mp.cpu_count(),inumsnaps)
	# 	nchunks=int(np.ceil(inumsnaps/float(nthreads)))
	# 	print("Using", nthreads,"threads to parse ",inumsnaps," snapshots in ",nchunks,"chunks")


	# 	for j in range(nchunks):
	# 		offset=j*nthreads + isnap
	# 		#if last chunk then must adjust nthreads
	# 		if (j==nchunks-1):
	# 			nthreads=inumsnaps-offset

	# 		snapKey = "Snap_%03d" %offset
	# 		ihaloset = np.where(MTFdata[snapKey]["HostHaloID"]>-1)[0] 

	# 		#when calling a process pass manager based proxies, which then are used to copy data back
	# 		processes=[mp.Process(target=SetHostIDParallel,args=(ihaloset,snap,MTFdata[snapKey]["HostHaloID"],MTFdata[snapKey]["origID"],mpMTFdata,HALOIDVAL)) for k in range(nthreads)]
	# 		#start each process


	# 		count=0
	# 		for p in processes:
	# 			p.start()
	# 			#space threads apart (join's time out is 0.25 seconds
	# 			p.join(0.2)
	# 			count+=1


	# 		# sel = np.where(MTFdata[snapKey]["HaloID"]>0)[0]
	# 		MTFdata[snapKey]["HostHaloID"][:] = np.array(mpMTFdata[snapKey][:])




	

	return MTFdata






				


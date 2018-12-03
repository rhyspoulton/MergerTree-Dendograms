import numpy as np 
import h5py
import convTreesToETF.VELOCIraptor_Python_Tools.velociraptor_python_tools as VPT

def convVELOCIraptorToMTF(opt,fieldsDict):

	treefields = ["ID","RootTail","Tail","Head","RootHead"]

	#Load in the VELOCIraptor catalogue
	Redshift,halodata,walkabletree = LoadVELOCIraptor(opt,treefields,fieldsDict.values())

	treedata = {}

	# Lets define a new dictionary for each snapshot
	for snap in range(opt.startSnap,opt.endSnap+1):

		snapKey = "Snap_%03d" %snap
		isnap = snap - opt.startSnap

		treedata[snapKey]  = {}

		#Loop over all the fields changing the keys to the MTF keys
		for field in fieldsDict.keys():

			if(fieldsDict[field] in treefields):
				treedata[snapKey][field] = walkabletree[isnap].pop(fieldsDict[field])

			#See if this dataset is the WWHalo Flag dataset
			elif(field=="WWHaloFlag"):
				tmpdata = halodata[isnap][fieldsDict[field]]

				#Create a boolean dataset the same shape as the data
				WWHaloFlag = np.zeros(tmpdata.shape,dtype=bool)

				# Mark where a WW halo is present
				WWHaloFlag[tmpdata==-1] = True

				#Delete the tmp data
				del tmpdata

				treedata[snapKey][field] = WWHaloFlag

			else:

				# Add the dataset into the treedata
				treedata[snapKey][field] = halodata[isnap].pop(fieldsDict[field])

	return Redshift,treedata


def LoadVELOCIraptor(opt,treefields,fieldKeys):

	numsnaps = opt.endSnap+1 - opt.startSnap

	extractfields = []
	for field in fieldKeys:

		#Skip if the field is in the treefields
		if(field in treefields):
			continue

		if(field=="Pos"):
			extractfields.extend(["Xc","Yc","Zc"])
		elif(field=="Vel"):
			extractfields.extend(["VXc","VYc","VZc"])
		else:
			extractfields.append(field)

	Redshift = np.zeros(numsnaps)

	halodata = [dict() for i in range(numsnaps)]

	#Read the VELOCIraptor properties files across the desired snapshots
	for snap in range(opt.startSnap,opt.endSnap+1):

		isnap = snap - opt.startSnap
		filename = opt.VELdir + "/snapshot_%03d.VELOCIraptor" %snap
		halodata[isnap],_ = VPT.ReadPropertyFile(filename,ibinary=2,desiredfields=extractfields)

		#Lets check if the position is in comoving units
		if(halodata[isnap]["UnitInfo"]["Comoving_or_Physical"]):

			#Convert to comoving 
			halodata[isnap]["Xc"]*=halodata[isnap]["SimulationInfo"]["h_val"]/halodata[isnap]["SimulationInfo"]["ScaleFactor"]
			halodata[isnap]["Yc"]*=halodata[isnap]["SimulationInfo"]["h_val"]/halodata[isnap]["SimulationInfo"]["ScaleFactor"]
			halodata[isnap]["Zc"]*=halodata[isnap]["SimulationInfo"]["h_val"]/halodata[isnap]["SimulationInfo"]["ScaleFactor"]

			#Lets convert all the types of radius
			for field in halodata[isnap].keys():
				if(field[0]=="R"):
					halodata[isnap][field]*=halodata[isnap]["SimulationInfo"]["h_val"]/halodata[isnap]["SimulationInfo"]["ScaleFactor"]

		Redshift[snap-opt.startSnap] = 1.0/halodata[isnap]["SimulationInfo"]["ScaleFactor"] - 1.0

		#Lets make sure all the units are in ETF

		#Distances in Mpc
		halodata[isnap]["Xc"]*=halodata[isnap]["UnitInfo"]["Length_unit_to_kpc"]/1000 #Mpc
		halodata[isnap]["Yc"]*=halodata[isnap]["UnitInfo"]["Length_unit_to_kpc"]/1000 #Mpc
		halodata[isnap]["Zc"]*=halodata[isnap]["UnitInfo"]["Length_unit_to_kpc"]/1000 #Mpc

		halodata[isnap]["Pos"] = np.column_stack([halodata[isnap].pop("Xc"),halodata[isnap].pop("Yc"),halodata[isnap].pop("Zc")])

		#Lets convert all the types of radius, velocity and masses
		for field in halodata[isnap].keys():
			if(field[0]=="R"):
				halodata[isnap][field]*=halodata[isnap]["UnitInfo"]["Length_unit_to_kpc"]/1000 #Mpc
			elif(field[0]=="M"):
				halodata[isnap][field]*=halodata[isnap]["UnitInfo"]["Mass_unit_to_solarmass"]/1e10 #1e10 solarmasses
			elif(field[0]=="V"):
				halodata[isnap][field]*=halodata[isnap]["UnitInfo"]["Velocity_unit_to_kms"] #1e10 solarmasses

		halodata[isnap]["Vel"] = np.column_stack([halodata[isnap].pop("VXc"),halodata[isnap].pop("VYc"),halodata[isnap].pop("VZc")])

	#Read in the walkable tree
	walkabletree,_ = VPT.ReadWalkableHDFTree(opt.VELwalkabletreefilename,False)



	return Redshift,halodata,walkabletree





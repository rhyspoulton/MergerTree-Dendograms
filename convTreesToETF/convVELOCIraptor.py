import numpy as np 
import h5py


def convVELOCIraptorToMTF(startSnap,endSnap,filename,fieldsDict):

	#Load in the VELOCIraptor catalogue
	Redshift,numsnaps,halodata = LoadVELOCIraptor(startSnap,endSnap,filename,fieldsDict.values())

	treedata = {}

	# Lets define a new dictionary for each snapshot
	for snap in range(startSnap,endSnap+1):

		snapKey = "Snap_%03d" %snap

		treedata[snapKey]  = {}

		#Loop over all the fields changing the keys to the MTF keys
		for field in fieldsDict.keys():

			#See if this dataset is the WWHalo Flag dataset
			if(field=="WWHaloFlag"):
				tmpdata = halodata[snapKey][fieldsDict[field]]

				#Create a boolean dataset the same shape as the data
				WWHaloFlag = np.zeros(tmpdata.shape,dtype=bool)

				# Mark where a WW halo is present
				WWHaloFlag[tmpdata==-1] = True

				#Delete the tmp data
				del tmpdata

				treedata[snapKey][field] = WWHaloFlag

			else:

				# Add the dataset into the treedata
				treedata[snapKey][field] = halodata[snapKey].pop(fieldsDict[field])

	return Redshift,treedata


def LoadVELOCIraptor(startSnap,endSnap,filename,fieldKeys):

	h=0.6751

	# Open up the  hdf file
	hdffile = h5py.File(filename,"r")

	# Extract the number of snapshots from the file
	numsnaps = endSnap-startSnap +1

	# Setup a dictionary to load the data into
	halodata = {}

	Redshift = np.zeros(numsnaps,dtype=float)

	# Read the header information
	icomove = hdffile["Header"]["Units"].attrs["Comoving_or_Physical"][...]
	length_unit_in_mpc = hdffile["Header"]["Units"].attrs["Length_unit_to_kpc"][...] / 1000.0
	mass_unit_in_10E10solarmass = hdffile["Header"]["Units"].attrs["Mass_unit_to_solarmass"][...] /1e10
	velocity_unit_in_kms = hdffile["Header"]["Units"].attrs["Velocity_unit_to_kms"][...]

	# Lets loop over snapshot and extract the necessary data sets
	for snap in range(startSnap,endSnap+1):

		snapKey = "Snap_%03d" %snap

		isnap = snap -startSnap

		Redshift[isnap] = 1.0/np.float(hdffile[snapKey].attrs["scalefactor"]) - 1.0

		# Every snapshot has a dictionary of the data
		halodata[snapKey] = {}

		#Loop over all the fields and extract them into memory
		for field in fieldKeys:
			
			#Do all the nessary conversions into ETF units
			if("Mass" in field):
				halodata[snapKey][field] = hdffile[snapKey][field][:] * mass_unit_in_10E10solarmass

			elif(field=="Pos"):

				if(icomove):
					halodata[snapKey][field] = np.column_stack([hdffile[snapKey]["Xc"][:],hdffile[snapKey]["Yc"][:],hdffile[snapKey]["Zc"][:]]) * length_unit_in_mpc
				else:
					halodata[snapKey][field] = np.column_stack([hdffile[snapKey]["Xc"][:],hdffile[snapKey]["Yc"][:],hdffile[snapKey]["Zc"][:]]) *h/np.float(hdffile[snapKey].attrs["scalefactor"]) * length_unit_in_mpc

			elif((field=="Radius") | (field[0:2]=="R_") | (field=="Rmax")) :

				if(icomove):
					halodata[snapKey][field] =hdffile[snapKey][field][:] * length_unit_in_mpc
				else:
					halodata[snapKey][field] =hdffile[snapKey][field][:] *h/np.float(hdffile[snapKey].attrs["scalefactor"]) * length_unit_in_mpc

			elif((field=="Vel") | (field=="Vmax") | (field=="sigV")):

				halodata[snapKey][field] = np.column_stack([hdffile[snapKey]["VXc"][:],hdffile[snapKey]["VYc"][:],hdffile[snapKey]["VZc"][:]]) * velocity_unit_in_kms

			else:
				halodata[snapKey][field] = hdffile[snapKey][field][:]



	hdffile.close()

	return Redshift,numsnaps,halodata





import sys
sys.path.append(sys.argv[0].replace("convToETF.py",""))
from convTreesToETF import *
from ETFio import WriteETFCatalogue
from ReadConfig import  ETFoptions
import h5py
import os

availableConverters = ["VEL","AHF","Rock","Dtrees"]

if(len(sys.argv)<3):
	raise SystemExit("No format of the merger tree parse. Usage: convToETF.py <merger tree format> <convToETF.cfg>")



if(sys.argv[1] not in availableConverters):
	raise SystemExit("Merger tree format %s not recognised, please choose from",availableConverters,"or create your own conversion tool to convert to ETF" %sys.argv[1])



opt = ETFoptions(sys.argv[1],sys.argv[2])


if(sys.argv[1]=="VEL"):
	####################################################################################


					#Converting VELOCIraptor & TreeFrog to ETF

	####################################################################################

	print("Doing the conversion of VELOCIraptor into ETF")

	fieldsDict = {}
	# The required Keys, do not change these!
	fieldsDict["HaloID"] = "ID"
	fieldsDict["StartProgenitor"] = "RootTail"
	fieldsDict["Progenitor"] = "Tail"
	fieldsDict["Descendant"] = "Head"
	fieldsDict["EndDescendant"] = "RootHead"
	fieldsDict["Pos"] = "Pos"
	fieldsDict["Vel"] = "Vel"
	fieldsDict["HostHaloID"] = "hostHaloID"
	fieldsDict["Mass"] = opt.MassDef
	fieldsDict["Radius"] = opt.RDef

	#Add in the extra fields into the fieldsDict
	for field in opt.ExtraFields:
		if field not in fieldsDict.values():
			fieldsDict[field] = field

	#Add an extra dataset so it is possible to tell if a halo is a WhereWolf halo
	if(opt.WWflag):
		fieldsDict["WWHaloFlag"] = "cNFW"

	print("Loading in the fields "+ " ".join(fieldsDict.values()))

	#The filename to read the data from
	Redshift,treedata = convVELOCIraptorToMTF(opt,fieldsDict)

elif(sys.argv[1]=="AHF"):
	####################################################################################


					#Converting AHF & MergerTree to MTF

	####################################################################################

	fieldsDict = {}
	# The required Keys, do not change these!
	fieldsDict["HaloID"] = ["ID(1)","int64"]
	fieldsDict["StartProgenitor"] = ["","int64"]
	fieldsDict["Progenitor"] = ["","int64"]
	fieldsDict["Descendant"] = ["","int64"]
	fieldsDict["EndDescendant"] = ["","int64"]
	fieldsDict["Pos"] = ["Xc(6),Yc(7),Zc(8)","float32"]
	fieldsDict["HostHaloID"] = ["hostHalo(2)","int64"]
	fieldsDict["Mass"] = [opt.MassDef,"float32"]
	fieldsDict["Radius"] = [opt.RDef,"float32"]

	# If the haloIDs aren't in the sussing format then the original ID's needs to be saved
	if(opt.sussingformat==False):
		fieldsDict["origID"] = ["ID(1)","int64"]

	if(len(opt.ExtraFields)!=len(opt.ExtraFieldsDtype)):
		raise SystemExit("Please input the datatype for each of the extra fields")

	#Add in the extra fields into the fieldsDict
	for field,dtype in zip(opt.ExtraFields,opt.ExtraFieldsDtype):
		if field not in fieldsDict:
			fieldsDict[field] = [field,dtype]

	print("Loading in the fields "+ " ".join([field[0] for field in fieldsDict.values()]))

	Redshift,treedata = convAHFToMTF(opt,fieldsDict)


elif(sys.argv[1]=="Rock"):
	####################################################################################


					#Converting Rockstar & consistent trees to MTF

	####################################################################################

	fieldsDict = {}
	# The required Keys, do not change these!
	fieldsDict["HaloID"] = ["id(1)","int64"]
	fieldsDict["StartProgenitor"] = ["","int64"]
	fieldsDict["Progenitor"] = ["","int64"]
	fieldsDict["Descendant"] = ["desc_id(3)","int64"]
	fieldsDict["EndDescendant"] = ["","int64"]
	fieldsDict["Pos"] = ["x(17),y(18),z(19)","float32"]
	fieldsDict["HostHaloID"] = ["upid(6)","int64"]
	fieldsDict["origID"] = ["id(1)","int64"]
	fieldsDict["Mass"] = [opt.MassDef,"float32"]
	fieldsDict["Radius"] = [opt.RDef,"float32"]

	if(len(opt.ExtraFields)!=len(opt.ExtraFieldsDtype)):
		raise SystemExit("Please input the datatype for each of the extra fields")

	#Add in the extra fields into the fieldsDict
	for field,dtype in zip(opt.ExtraFields,opt.ExtraFieldsDtype):
		if field not in fieldsDict:
			fieldsDict[field] = [field,dtype]

	print("Loading in the fields "+ " ".join([field[0] for field in fieldsDict.values()]))

	# Extra fields required to build the ETF catalogue
	snapColName = "Snap_idx(31)"
	numProgName = "num_prog(4)"
	redshiftColName= "scale(0)"


	Redshift,treedata = convRockstarToMTF(snapColName,numProgName,redshiftColName,opt.startSnap ,opt.endSnap,opt.Rockfilelist,fieldsDict)


elif(sys.argv[1]=="Dtrees"):
	####################################################################################


					#Converting D-trees format

	####################################################################################

	fieldsDict = {}
	# The required Keys, do not change these!
	fieldsDict["HaloID"] = "/haloTrees/nodeIndex"
	fieldsDict["Descendant"] = "/haloTrees/descendantIndex"
	fieldsDict["DescendantSnap"] = "/haloTrees/descendantSnapshot"
	fieldsDict["Progenitor"] = "/haloTrees/mainProgenitorIndex"
	fieldsDict["Mass"] = opt.MassDef
	fieldsDict["Radius"] = opt.RDef
	fieldsDict["Pos"] = "/haloTrees/position"
	fieldsDict["Vel"] = "/haloTrees/velocity"
	fieldsDict["AngMom"] = "/haloTrees/angularMomentum"
	fieldsDict["HostHaloID"] = "/haloTrees/hostIndex"
	fieldsDict["isMainProgenitor"] = "/haloTrees/isMainProgenitor"

	#Add in the extra fields into the fieldsDict
	for field in opt.ExtraFields:
		fieldsDict[field] = field

	# Extra fields required to build the ETF catalogue
	snapKey = "/haloTrees/snapshotNumber"
	scalefactorKey = "/outputTimes/redshift"

	treedata,Redshift = convDtreesToMTF(opt,snapKey,scalefactorKey,fieldsDict)

	fieldsDict["StartProgenitor"] = ""
	fieldsDict["EndDescendant"] = ""
	fieldsDict["DescendantIndex"] = ""
	fieldsDict["DescendantSnap"] = ""
	fieldsDict["ProgenitorIndex"] = ""
	fieldsDict["ProgenitorSnap"] = ""


WriteETFCatalogue(sys.argv[1],opt,treedata,Redshift,fieldsDict)





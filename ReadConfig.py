import numpy as np



class plotOptions(object):

	def __init__(self,filename,outdir):

		#Lets set the defaults
		self.outdir = outdir
		self.maxNumBranches = 20
		self.maxdepth = 2
		self.plotSubhaloBranches = 1
		self.minNsnapsExist = 10
		self.showBranchTypeLabel = 1
		self.sizeLabel = "M$_{\rm vir}$"
		self.sizeUnit = "10$^{10}$ M$_{\odot}$"
		self.maxSizeFontSize = 24
		self.maxSizeFormat = "%.2f"
		self.maxSizeFontDist = 10
		self.logged = 0
		self.subBranchSizeFactor = 1.0
		self.plotNumRvir = 2.5
		self.xLabel = "R$_{\rm#orbit}$/R$_{\rm#vir,parent}$"
		self.marker = "line"
		self.numSubplotsMain = 4
		self.plotColorBar = 0
		self.colMap = "winter"
		self.cbarLabel="R$_{\rm vir}$#[Mpc]"
		self.fileDesc = "dend"
		self.snapoffset = 0
		self.overplotdata = 0
		self.insetPlot = 1
		self.WWflag = 0

		with open(filename,"r") as f:

			for line in f:

				if(line[0]=="#"): continue

				line = line.replace(" ","")

				line = line.strip()

				if(not line): continue

				line = line.split("=")


				if(line[0]=="maxNumBranches"):
					self.maxNumBranches = int(line[1])

				elif(line[0]=="maxdepth"):
					self.maxdepth = int(line[1])

				elif(line[0]=="plotSubhaloBranches"):
					self.plotSubhaloBranches = int(line[1])

				elif(line[0]=="minNsnapsExist"):
					self.minNsnapsExist = int(line[1])

				elif(line[0]=="showBranchTypeLabel"):
					self.showBranchTypeLabel = int(line[1])

				elif(line[0]=="sizeLabel"):
					self.sizeLabel = line[1].replace("#"," ")

				elif(line[0]=="sizeUnit"):
					self.sizeUnit = line[1].replace("#"," ")

				elif(line[0]=="maxSizeFormat"):
					self.maxSizeFormat = line[1]

				elif(line[0]=="maxSizeFontSize"):
					self.maxSizeFontSize = int(line[1])

				elif(line[0]=="maxSizeFontDist"):
					self.maxSizeFontDist = int(line[1])

				elif(line[0]=="logged"):
					self.logged = int(line[1])

				elif(line[0]=="subBranchSizeFactor"):
					self.subBranchSizeFactor = float(line[1])

				elif(line[0]=="plotNumRvir"):
					self.plotNumRvir = float(line[1])

				elif(line[0]=="xLabel"):
					self.xLabel = line[1].replace("#"," ")

				elif(line[0]=="marker"):
					self.marker = line[1]

				elif(line[0]=="numSubplotsMain"):
					self.numSubplotsMain = int(line[1])

				elif(line[0]=="plotColorBar"):
					self.plotColorBar = int(line[1])

				elif(line[0]=="colMap"):
					self.colMap = line[1]

				elif(line[0]=="cbarLabel"):
					self.cbarLabel = line[1].replace("#"," ")

				elif(line[0]=="fileDesc"):
					self.fileDesc = line[1]

				elif(line[0]=="snapoffset"):
					self.snapoffset = int(line[1])

				elif(line[0]=="overplotdata"):
					self.overplotdata = int(line[1])

				elif(line[0]=="overplotFormat"):
					self.overplotFormat = line[1]

				elif(line[0]=="insetPlot"):
					self.insetPlot = int(line[1])

				elif(line[0]=="WWflag"):
					self.WWflag = int(line[1])

				else:
					raise OSError("Invalid config option %s, please only use the options in the sample config file" %line[0])


class ETFoptions(object):

	def __init__(self,MTF,filename):

		#Set the defaults
		self.startSnap = 0
		self.endSnap = 0
		self.Nsnaps = 0
		self.h = 0.6751
		self.boxsize = 0
		self.HALOIDVAL = 1000000000000
		self.outfilename = MTF
		self.MassDef = ""
		self.RDef = ""
		self.ExtraFields = []
		self.ExtraFieldsDtype = []
		self.VELdir = ""
		self.VELwalkabletreefilename = ""
		self.WWflag = 0
		self.AHFhalofilelist = ""
		self.AHFtreefilelist = ""
		self.sussingformat = 0
		self.Rockfilelist = ""
		self.Millfilename = ""
		self.iverbose = 0

		with open(filename,"r") as f:

			for line in f:

				if(line[0]=="#"): continue

				line = line.replace(" ","")

				line = line.strip()

				if(not line): continue

				line = line.split("=")

				if(line[0]=="startSnap"):
					self.startSnap = int(line[1])

				elif(line[0]=="endSnap"):
					self.endSnap = int(line[1])

				elif(line[0]=="Nsnaps"):
					self.Nsnaps = int(line[1])

				elif(line[0]=="h"):
					self.h = float(line[1])

				elif(line[0]=="boxsize"):
					self.boxsize = float(line[1])

				elif(line[0]=="HALOIDVAL"):
					self.HALOIDVAL = np.uint64(line[1])

				elif(line[0]=="outfilename"):
					self.outfilename = line[1]

				elif(line[0]=="MassDef"):
					self.MassDef = line[1]

				elif(line[0]=="RDef"):
					self.RDef = line[1]

				elif(line[0]=="ExtraFields"):
					if(line[1]):
						self.ExtraFields = line[1].split(",")


				elif(line[0]=="ExtraFieldsDtype"):
					if(line[1]):
						self.ExtraFieldsDtype = line[1].split(",")

				# VELOCIraptor specifics

				elif(line[0]=="VELdir"):
					if(MTF=="VEL"):
						self.VELdir = line[1]

				elif(line[0]=="VELwalkabletreefilename"):
					if(MTF=="VEL"):
						self.VELwalkabletreefilename = line[1]

				elif(line[0]=="WWflag"):
					self.WWflag = int(line[1])

				# AHF specifics

				elif(line[0]=="AHFhalofilelist"):
					if(MTF=="AHF"):
						self.AHFhalofilelist = line[1]

				elif(line[0]=="AHFtreefilelist"):
					if(MTF=="AHF"):
						self.AHFtreefilelist = line[1]

				elif(line[0]=="sussingformat"):
					if(MTF=="AHF"):
						self.sussingformat = int(line[1])

				# Rockstar specifics

				elif(line[0]=="Rockfilelist"):
					if(MTF=="Rock"):
						self.Rockfilelist = line[1]

				# Millenium specifics

				elif(line[0]=="Millfilename"):
					if(MTF=="Mill"):
						self.Millfilename = line[1]

				elif(line[0]=="iverbose"):
					self.iverbose=int(line[1])

				else:
					raise OSError("Invalid config option %s, please only use the options in the sample config file" %line[0])

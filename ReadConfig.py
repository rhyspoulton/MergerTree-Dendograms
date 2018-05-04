import numpy as np



class plotOptions(object):

	def __init__(self,filename,outdir):

		self.outdir = outdir

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

				elif(line[0]=="sizePoint"):
					self.sizePoint = float(line[1])

				elif(line[0]=="sizeLabel"):
					self.sizeLabel = line[1].replace("#"," ")

				elif(line[0]=="maxSizeFormat"):
					self.maxSizeFormat = line[1]

				elif(line[0]=="logged"):
					self.logged = int(line[1])

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

				else:
					raise OSError("Invalid config option %s, please only use the options in the sample config file" %line[0])


class ETFoptions(object):

	def __init__(self,MTF,filename):

		self.headerDict = {}

		with open(filename,"r") as f:

			for line in f:

				if(line[0]=="#"): continue

				line = line.replace(" ","")

				line = line.strip()

				if(not line): continue

				line = line.split("=")

				if(line[0]=="startSnap"):
					self.headerDict["startSnap"] = int(line[1])

				elif(line[0]=="endSnap"):
					self.headerDict["endSnap"] = int(line[1])

				elif(line[0]=="Nsnaps"):
					self.headerDict["Nsnaps"] = int(line[1])

				elif(line[0]=="Munit"):
					self.headerDict["Munit"] = float(line[1])

				elif(line[0]=="h"):
					self.headerDict["h"] = float(line[1])

				elif(line[0]=="boxsize"):
					self.headerDict["boxsize"] = float(line[1])

				elif(line[0]=="HALOIDVAL"):
					self.headerDict["HALOIDVAL"] = np.uint64(line[1])

				elif(line[0]=="outfilename"):
					self.outfilename = line[1]

				elif(line[0]=="MassDef"):
					self.MassDef = line[1]

				elif(line[0]=="RDef"):
					self.RDef = line[1]

				elif(line[0]=="ExtraFields"):
					if(line[1]):
						self.ExtraFields = line[1].split(",")
					else:
						self.ExtraFields = []

				elif(line[0]=="ExtraFieldsDtype"):
					if(line[1]):
						self.ExtraFieldsDtype = line[1].split(",")
					else:
						self.ExtraFieldsDtype = []

				# VELOCIraptor specifics

				elif(line[0]=="VELfilename"):
					if(MTF=="VEL"):
						self.VELfilename = line[1]

				# AHF specifics

				elif(line[0]=="AHFhalofilelist"):
					if(MTF=="AHF"):
						self.AHFhalofilelist = line[1]

				elif(line[0]=="AHFtreefilelist"):
					if(MTF=="AHF"):
						self.AHFtreefilelist = line[1]

				# Rockstar specifics

				elif(line[0]=="Rockfilelist"):
					if(MTF=="Rock"):
						self.Rockfilelist = line[1]

				# Millenium specifics

				elif(line[0]=="Millfilename"):
					if(MTF=="Mill"):
						self.Millfilename = line[1]

				else:
					raise OSError("Invalid config option %s, please only use the options in the sample config file" %line[0])

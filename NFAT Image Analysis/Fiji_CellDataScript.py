import sys
import os
from ij import IJ
from ij.measure import ResultsTable
from ij.gui import WaitForUserDialog
from ij import WindowManager
from ij.io import DirectoryChooser
from ij.gui import GenericDialog
from ij.gui import Roi
from ij.gui import Overlay
from ij.gui import TextRoi
from ij.plugin.frame import RoiManager

from ij.gui import WaitForUserDialog
from ij.process import AutoThresholder
from fiji.threshold import Auto_Threshold

from java.io import File
from java.awt import Color

reload(sys)
sys.setdefaultencoding('utf-8')

IJ.log("test")
IJ.run("Set Measurements...", "area mean perimeter integrated redirect=None decimal=3");
#DC = DirectoryChooser("choose trackmate data export directory")
#path = DC.getDirectory()

"""
USER INPUT SEGMENT
Please personalize Value here:

note for batch processing: please organize each dataset in a folder and put all of this folder inside another one.
Choose this folder to batch process all dataset folders inside.

If False - please select only one folder containing tifs for processing
"""
batch = False

autoSaveFile = True								
exportFolder = "Path/To/Export/Folder/Results"				#Where to save	

#Auto threshold Methods. Copy from macro recorder	
thMethodeCell = "Li dark"															
thMethodeDapi = "Li dark"		
thMethodeNFAT = "Li dark"		

#settings for image preprocessing before autothreshold
dapiFilter = "Median..."    		#uncomment which you want to use						
filter = "Gaussian Blur..."
uSigma = str(2)					#sigma for filter, pixels

#Channels
actinChannel = 1
NFATChannel = 2
DAPIChannel = 3

#Background
doRollingBall = True
rollingBallRadius = 20 #px

datashort = "ome"
filt = "companion"

#Z Slice
zslice = 7

#Export Proof Images with created Rois and example Data for each cell
debug = True
debugFolderName = "debug"	#Results subfolder

#Options for creation of cytoplasma Rois. 
saveOuterRois = False			#set to true for creation of Roi sets for new datasets. Not recommended in batch mode
loadOuterRois = True			#only possible when Rois have already been created. Macro will break else. 
roiFolderName = "cytoRoi"

"""
END OF USER SEGMENT
"""


IJ.setForegroundColor(255, 255, 0)
if debug:
	if not debugFolderName in os.listdir(exportFolder):
		os.mkdir(exportFolder + "/" + debugFolderName)
		


aT = AutoThresholder()


dc = DirectoryChooser("choose directory")

directory = " " 
directory = dc.getDirectory()

if batch:
	folders = os.listdir(directory)
else: 
	folders = [""]
for folder in folders:

	files = os.listdir(directory + folder)
	files = [file for file in files if filt in file]
	
	image = 0
	Results = ResultsTable()
	for file in files:
		
		
		rM = RoiManager.getRoiManager()
		rM.reset()
		
		IJ.log(file)
		#imp = IJ.openImage(directory+folder+file);
		#imp.show()
		IJ.run("Bio-Formats Importer", "open=["+directory+folder + "/" +file+"]" + " color_mode=Grayscale display_rois rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		imp = IJ.getImage()
		
		 
		if debug:
			dbImp = imp.duplicate()
			dbIp = dbImp.getProcessor()
			dbOverlay = Overlay()

	
		IJ.log("checkpoint 1")
		imp.setZ(zslice)
		imp.setC(actinChannel)
		IJ.run(imp, "Enhance Contrast", "saturated=0.35");
		IJ.run(imp, filter , "sigma=" + uSigma + " slice");
		
		
		gotRois = True
		if loadOuterRois:
			roiFile = file.replace("." + datashort, ".zip")
			if roiFile in os.listdir(directory + folder + "/" + roiFolderName):
				rM.open(directory + folder + "/" + roiFolderName + "/" + roiFile)
			else: gotRois = False
		else:
			IJ.setAutoThreshold(imp, thMethodeCell);
			IJ.run(imp, "Analyze Particles...", "size=10-Infinity display exclude clear include add slice")
			WJ = WaitForUserDialog("CheckRois")
			WJ.show()
		
		
		if saveOuterRois:
			if roiFolderName in os.listdir(directory+folder):
				roiFile = file.replace("." + datashort, ".zip")
				rM.save(directory + folder + "/" + roiFolderName + "/" +roiFile)
			else: 
				os.mkdir(directory + folder +"/" + roiFolderName)
				roiFile = file.replace("." + datashort, ".zip")
				rM.save(directory + folder + "/" + roiFolderName + "/" +roiFile)
				
		if gotRois:	
			
			Rois = rM.getRoisAsArray()
			IJ.log("Rois Array length: " + str(len(Rois)))
			roiInd = 0
			
			dapiRois = []
			
			
			imp.setC(NFATChannel)
			
			#histogram = imp.getProcessor().convertToByte(True).getStats().histogram()
			#histogram = [int(dat) for dat in histogram]
			#threshold = aT.getThreshold("Li", histogram)
			
			
			imp.setC(DAPIChannel) 
			IJ.run(imp, dapiFilter , "sigma=" + uSigma + " slice");
			IJ.setAutoThreshold(imp, thMethodeDapi);
			roiSets = []
			
			for roi in Rois:
					
				IJ.log("DapiRoi")
				imp.setRoi(roi)
				rM.reset()
				
				IJ.run(imp, "Analyze Particles...", "size=10-Infinity display exclude clear include add slice");
				dapiRoisTemp = rM.getRoisAsArray()
				
				if dapiRoisTemp:
					
					dapiRoi = dapiRoisTemp[0]
					roiSets.append([roi,dapiRoi])
				
			IJ.log("DapiQuit")
			IJ.log(str(roiSets))
			roiInd = 0
			for roiSet in roiSets:
				IJ.log("roiSet")
				roi = roiSet[0]
				dapiRoi = roiSet[1]
				imp.setC(NFATChannel)
				imp.setRoi(roi)
				overall_int_dens = IJ.getValue(imp, "RawIntDen")
				
				imp.setRoi(dapiRoi)
				nucleus_int_dens = IJ.getValue(imp, "RawIntDen")
				nucleus_mean = IJ.getValue(imp, "Mean")
				
				cytopRoi = Roi.xor([roi, dapiRoi])
				imp.setRoi(cytopRoi)
				cytoplasma_mean = IJ.getValue(imp, "Mean")
				
				Results.addRow()
				Results.addValue("Sample", file.split(".")[0])
				Results.addValue("Image", image)
				Results.addValue("Cell Tag", roiInd)
				Results.addValue("overall integrated density", overall_int_dens)
				Results.addValue("nucleus integrated density", nucleus_int_dens)
				Results.addValue("cytoplasma mean", cytoplasma_mean)
				Results.addValue("nucleus mean", nucleus_mean)
				
				if debug:
				
					IJ.log("test")
					bounds = roi.getBounds()
					nfatRatio = nucleus_int_dens/overall_int_dens
					bleftx = bounds.x
					blefty = bounds.y + bounds.height
					#dp.drawString(str(nfatRatio), bleftx, blefty);
					
					dbOverlay.add(TextRoi(bleftx, blefty, "Cell: " +  str(roiInd) + "\n NFAT Ratio: " + str(nfatRatio)))
					
				roiInd = roiInd + 1
				
				
				
				
			imp.changes = False
			imp.close()
			
			if debug:
		
				dbOverlay.add(TextRoi(10, 10, "Settings: " + "\n"
				"Threshhold Methode Dapi" + thMethodeDapi + "\n" +
				"filter: " + filter + "\n" + 
				"sigma: " + str(uSigma)))
				dbImp.setC(actinChannel)
				IJ.run(dbImp, "Enhance Contrast", "saturated=0.35");
				dbImp.setC(NFATChannel)
				IJ.run(dbImp, "Enhance Contrast", "saturated=0.35");
				dbImp.setC(DAPIChannel)
				IJ.run(dbImp, "Enhance Contrast", "saturated=0.35");
				dbImp.setZ(zslice)
				if not batch:
					dbImp.show()
				dbImp.setOverlay(dbOverlay)
				
				for roiSet in roiSets:
					roi = roiSet[0]
					dapiRoi = roiSet[1]
	
					roi.setPosition(0)
					roi.setStrokeColor(Color.yellow)
					roi.setStrokeWidth(1.0)
					roi.setStrokeWidth(1.0)
					dbOverlay.add(roi)
	
					dapiRoi.setPosition(0)
					dapiRoi.setStrokeColor(Color.magenta)
					dapiRoi.setStrokeWidth(1.0)
					dbOverlay.add(dapiRoi)
				
				
				IJ.saveAs(dbImp, "Tiff", exportFolder + "/debug/" + file);
				
				if batch:
					dbImp.changes = False
					dbImp.close()
				
			image = image+1
			
		else: 
			IJ.log("No Rois")
			imp.changes = False
			imp.close()
	
	
	Results.show("End Results")
	if autoSaveFile:
		
		filename = file.split(".")[0] + ".csv"
		Results.save(exportFolder + "/" + filename)
	
	
	
	



	
	
	
	
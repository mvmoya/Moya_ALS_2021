## Tile_scan_stitching v0.1
## Developed by VMoya 06/02/2018 (last update: 06/08/2018)
## Welcome to Tile_scan_stitching!
## This program works best when run on a RGB composite image that was acquired as a tile scan
## on the LSM 700 in the Heintz Lab. Pixel height data is key for the default seam overlays
## to be representative of the likely seams, but the program offers the opportunity to set the
## width of the seams manually should this automated part fail. 100X images are not currently
## supported.
## Again, the program works best when pixel resolution was carried over in microns. If not,
## the estimation of the seams will be totally off.
## Good luck and happy stitching!
## This is an ImageJ plugin

import ij.plugin as plugin
import ij.gui as gui
import ij.process
from ij import ImageStack
from ij import IJ
import ij.io.FileInfo as fileInfo

imp = IJ.getImage()
impProcessor = IJ.getProcessor()
impProcessor.setLineWidth(1)
impProcessor.setColor(25)
imp.setOverlay(gui.Overlay())
pixelHeight = round(imp.getLocalCalibration().pixelHeight, 4)
name = imp.getShortTitle()
upperName = name.upper()
magnification = ''


if [x in upperName for x in ['5X', '10X', '20X', '40X', '63X']].count(True) == 1:

	window = gui.NonBlockingGenericDialog("Specify Tile Scan Properties")
	window.addMessage("What are the tile scan dimensions?")
	window.addNumericField('Height:  ', 2.0, 1)
	window.addNumericField('Width:  ', 3.0, 1)
	window.addCheckbox('Set seam sizes manually?   ', False)
	window.addCheckbox('Approve seams first?   ', True)
	window.showDialog()
	
	if window.wasCanceled():
	
		raise Exception("Try Again")
		
	yTiles = window.getNextNumber()
	xTiles = window.getNextNumber()
	manually = window.getNextBoolean()
	approve = window.getNextBoolean()

else:

	window = gui.NonBlockingGenericDialog('Specify Tile Scan Properties')
	window.addMessage('What are the tile scan dimensions?')
	window.addNumericField('Height:  ', 2.0, 1)
	window.addNumericField('Width:  ', 3.0, 1)
	window.addChoice('Magnification:  ', ['5X', '10X', '20X', '40X', '63X'], '10X')
	window.addCheckbox('Set seam sizes manually?   ', False)
	window.addCheckbox('Approve seams first?   ', True)
	window.showDialog()
	
	if window.wasCanceled():
	
		raise Exception('Try Again')
		
	yTiles = window.getNextNumber()
	xTiles = window.getNextNumber()
	magnification = window.getNextChoice()
	manually = window.getNextBoolean()
	approve = window.getNextBoolean()
	
if xTiles != xTiles or yTiles != yTiles:

	errorDialog = gui.NonBlockingGenericDialog('Error')
	errorDialog.addMessage('Dimensions should be expressed as valid integers')
	errorDialog.hideCancelButton()
	errorDialog.showDialog()
	raise ValueError('Dimensions should be expressed as valid integers')

imageHeight = imp.getFileInfo().height
tileHeight = imageHeight / yTiles
imageWidth = imp.getFileInfo().width
tileWidth = imageWidth / xTiles
seamOverlay = gui.Overlay()


if not manually:

	factor = pixelHeight / 0.3126 # This seems to be the smallest pixel height in microns for our standard tile scans? Lies.
	
	if '10X' in upperName or magnification == '10X':

		seamHeight = 14 / factor
		seamWidth = 10 / factor
		
	elif '20X' in upperName or magnification == '20X':

		seamHeight = 12 / factor
		seamWidth = 0 / factor
	
	elif '40X' in upperName or magnification == '40X':
	
		seamHeight = 8 / factor
		seamWidth = 6 / factor
		
	elif '63X' in upperName or magnification == '63X':
		seamHeight = 3 / factor
		seamWidth = 0 / factor
	
	elif '5X' in upperName or magnification == '5X':
	
		seamHeight = 24 / factor
		seamWidth = 0 / factor
	
	else:

		errorDialog = gui.NonBlockingGenericDialog('Oops, Sorry!')
		errorDialog.addMessage("She hasn't added this functionality yet!")
		errorDialog.hideCancelButton()
		errorDialog.showDialog()
		raise Exception('Check back later! GTFO')

	for i in range(int(yTiles) - 1):

		hSeam = gui.Roi(0, tileHeight + tileHeight * i, imageWidth, seamHeight)
		seamOverlay.add(hSeam)

		for j in range(int(xTiles) - 1):
	
			vSeam = gui.Roi(tileWidth + tileWidth * j, 0, seamWidth, imageHeight)
			seamOverlay.add(vSeam)

	imp.setOverlay(seamOverlay)
	window = gui.NonBlockingGenericDialog("Approve the seams")
	window.addMessage("Zoom in and check that the seams\r\nare in the correct places.\r\nThen click OK.\r\nIf wrong, click Adjust.")
	window.setCancelLabel('Adjust')
	window.showDialog()
	
	if window.wasCanceled():
	
		manually = True

while manually:

	try:
		seamHeight = float(seamHeight)
		seamWidth = float(seamWidth)
		
	except:
		seamHeight = 7.0
		seamWidth = 0.0
	
	window = gui.NonBlockingGenericDialog('Set Seam Sizes Manually')
	window.addMessage('Set the sizes of the seams in pixels')
	window.addNumericField('Horizontal seam (px):   ', seamHeight, 1)
	window.addNumericField('Vertical seam (px):   ', seamWidth, 1)
	window.setOKLabel('Update seams')
	window.showDialog()
	
	if window.wasCanceled():
	
		raise Exception('Was Canceled')
		break
		
	seamHeight = window.getNextNumber()
	seamWidth = window.getNextNumber()
	
	if seamHeight != seamHeight or seamWidth != seamWidth:
	
		errorDialog = gui.NonBlockingGenericDialog()
		errorDialog.addMessage('Seam sizes must be valid integers')
		errorDialog.hideCancelButton()
		errorDialog.showDialog()
		raise ValueError('Seam sizes must be valid integers')
		break

	imp.setOverlay(gui.Overlay())
	seamOverlay = gui.Overlay()
	
	for i in range(int(yTiles) - 1):

		hSeam = gui.Roi(0, tileHeight + tileHeight * i, imageWidth, seamHeight)
		seamOverlay.add(hSeam)

		for j in range(int(xTiles) - 1):
	
			vSeam = gui.Roi(tileWidth + tileWidth * j, 0, seamWidth, imageHeight)
			seamOverlay.add(vSeam)

	imp.setOverlay(seamOverlay)
	window = gui.NonBlockingGenericDialog("Approve the seams")
	window.addMessage("Zoom in and check that the seams\r\nare in the correct places.\r\nThen click OK.\r\nIf wrong, click Cancel.")
	window.setCancelLabel('Adjust')
	window.showDialog()
	
	if window.wasCanceled():
	
		continue

	elif window.wasOKed():
		
		break

for i in range(int(yTiles)):

	if i == 0:	
		yMin = seamHeight
		
	else:	
		yMin = (tileHeight * i) + seamHeight
	
	yMax = tileHeight - seamHeight
	
	for j in range(int(xTiles)):	
		
		if j == 0:		
			xMin = seamWidth
			
		else:		
			xMin = (tileWidth * j) + seamWidth
		
		xMax = tileWidth - seamWidth

		if i == 0 and j == 0:		
			finalStack = ImageStack(int(xMax), int(yMax))
		
		tempRoi = gui.Roi(xMin, yMin, xMax, yMax)
		roiImage = imp.setRoi(tempRoi, False)
		duped = imp.duplicate()
		dupedProcessor = duped.getProcessor()
		finalStack.addSlice(dupedProcessor)

stackImage = ij.ImagePlus('Stack', finalStack)
montageImage = plugin.MontageMaker().makeMontage2(stackImage, int(xTiles), int(yTiles), 1.0, 1, int(xTiles*yTiles), 1, 0, False)
montageImage.setTitle('Stitched_'+name)
montageImage.show()

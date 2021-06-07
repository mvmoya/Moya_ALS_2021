## This program calculates the distance of ROI points in an image to an ROI line (or two).
## It does this by generating a perpendicular line that can be made to pass through the ROI point.
## It then calculates the distance between your point and the intersection of this perpendicular
## line with the reference ROI line. The functions called in the program are at the top.
## Please note that before running the program, the first ROI line must already be drawn. Have fun.
## Java-Python (Jython) sucks.

import java.lang.Math as math
import ij.measure as measure
import ij.plugin
from java.awt import Color
from java.awt import Frame
from java.awt import Font
from java.lang import String
from java.awt.event import MouseAdapter
from java.awt.event import TextListener
from ij import IJ
import ij.gui as gui

def Get_Line_Info(line):
	angle = line.getAngle()
	xbase = line.getXBase()
	ybase = -line.getYBase()
	width = line.getFloatWidth()
	height = line.getFloatHeight()
			
	if angle < 0:
		xend = xbase + width
		yend = ybase - height
		slope = -(height/width)
	else:	
		ybase = ybase - height
		xend = xbase + width
		yend = -line.getYBase()
		slope = (height/width)

	return (angle, xbase, ybase, width, height, xend, yend, slope)

def Perpendicular_line(imp, line_info):
	angle, xbase, ybase, width, height, xend, yend, slope = line_info
	
	width, height, nChannels, nSlices, nFrames = imp.getDimensions()
	b = ybase - (slope*xbase)
	perp_angle = math.toRadians(angle - 90)
	slope_new_line = math.tan(perp_angle)
			
	return slope_new_line, b

def Distance_point_to_line(line_info, new_slope, b, cell_points):
	angle, xbase, ybase, width, height, xend, yend, slope = line_info
	cell_distance = []
	for index in range(len(cell_points)):
		x_coor, y_coor = cell_points[index]
		y_coor = -y_coor
		b_new_line = y_coor - (new_slope*x_coor)
		x_intercept = (b - b_new_line)/(new_slope - slope)
		y_intercept = (new_slope*x_intercept + b_new_line)
		#print(b_new_line, x_intercept, y_intercept)
		distance = ((x_intercept- x_coor)**2 + (y_intercept - y_coor)**2)**0.5
		cell_distance.append(distance)

	return cell_distance

def Generate_table(table, cell_distance_pia, cell_distance_wm, cell_points):
	
	for i in range(len(cell_distance_pia)):
		xcoor, ycoor = cell_points[i]
		table.incrementCounter()
		table.addValue("Cell number", "Cell # " + str(i + 1))
		table.addValue("X coordinate", xcoor)
		table.addValue("Y coordinate", ycoor)
		table.addValue("Distance from Pia", cell_distance_pia[i])
		table.addValue("Distance from WM", cell_distance_wm[i])
		pia = cell_distance_pia[i]
		wm = cell_distance_wm[i]
		percent_distance = (pia/(wm + pia))*100
		table_data.addValue("Percent depth from Pia", percent_distance)
				
	table_data.show("Results in Pixels")

table_data = measure.ResultsTable()

imp = IJ.getImage()
line_pia = imp.getRoi()
overlay = gui.Overlay(line_pia)

if line_pia is None:
	window = gui.NonBlockingGenericDialog("Pial Surface")
	window.addMessage("Please draw straight line to represent pia.\r\nThen click OK")
	window.showDialog()
	if window.wasCanceled():
		raise NameError("Try again")

line_info_pia = Get_Line_Info(line_pia)

window = gui.NonBlockingGenericDialog("Draw white matter line")
window.addMessage("Draw the next line to outline white matter.\r\nThen click OK")
window.showDialog()
if window.wasCanceled():
	raise NameError("Try again")

line_wm = imp.getRoi()
overlay.add(line_wm)
line_info_wm = Get_Line_Info(line_wm)

window = gui.NonBlockingGenericDialog("Cell selection")
window.addMessage("Using the multi-point tool, select each cell to be analyzed.\r\nThen click OK")
window.showDialog()
if window.wasCanceled():
	raise NameError("Try again")

cells = imp.getRoi()
overlay.add(cells)
xpoints = cells.getXCoordinates()
ypoints = cells.getYCoordinates()

xbase = cells.getXBase()
ybase = cells.getYBase()

cell_coordinates = []
for item in xpoints:
	n = xpoints.index(item)
	x = xpoints[n] + xbase
	y = ypoints[n] + ybase
	coordinate = (x, y)
	cell_coordinates.append(coordinate)
	


new_slope_pia, b = Perpendicular_line(imp, line_info_pia)
cell_distance_from_pia = Distance_point_to_line(line_info_pia, new_slope_pia, b, cell_coordinates)
new_slope_wm, b = Perpendicular_line(imp, line_info_wm)
cell_distance_from_wm = Distance_point_to_line(line_info_wm, new_slope_pia, b, cell_coordinates)
Generate_table(table_data, cell_distance_from_pia, cell_distance_from_wm, cell_coordinates)

imp.setOverlay(overlay)

new_rgb = imp.flatten()
new_rgb.show()

#print(cell_coordinates)

#print(new_slope_pia)
#print(new_slope_wm)
#print(b)
#print(cell_distance_from_pia)
#print(cell_distance_from_wm)


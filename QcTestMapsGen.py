import arcpy, os
from arcpy import env  
from arcpy.sa import *

# Set environment settings
dataFrame = arcpy.mapping.MapDocument("CURRENT").activeDataFrame			## Define the dataFrame


## STEP 1: CREATION OF THE TABLE ---------------------------------------- 
table = arcpy.GetParameter(0)									## Load the table
tabName = arcpy.GetParameterAsText(1)							## Set the table name  Tik_1km
outFolder = "C:/Roby/PhD/5. Deception/GIS/" + tabName			## Destination folder  "C:/Roby/PhD/5. Deception/GIS/Tik_1km"
filePath = outFolder + "/" + tabName + ".shp"					## Path to file "C:/Roby/PhD/5. Deception/GIS/Tik_1km.shp"

# Set the local variables
in_Table = table		
x_coords = "Lon_x"
y_coords = "Lat_y"
z_coords = "Dep_z"
out_Layer = tabName

# Set the spatial reference
spRef = r"Coordinate Systems\Projected Coordinate Systems\Utm\WGS 1984\Southern Hemisphere\WGS 1984 UTM Zone 20S.prj"
 
# Make the XY event layer...
arcpy.MakeXYEventLayer_management(in_Table, x_coords, y_coords, out_Layer, spRef, z_coords)
 
# Execute FeatureClassToGeodatabase
if (os.path.isdir(outFolder)):
    arcpy.FeatureClassToShapefile_conversion(out_Layer, outFolder)
else:
    arcpy.CreateFolder_management("C:/Roby/PhD/5. Deception/GIS/", tabName) ## Create destination folder
	
arcpy.FeatureClassToShapefile_conversion(out_Layer, outFolder)

## END STEP 1 ----------------------------------------

tabLayer = arcpy.mapping.Layer(filePath)
arcpy.mapping.AddLayer(dataFrame,tabLayer,"BOTTOM")					## Add the layer to the active dataFrame Tik_1km.shp

## STEP 2: MAPS CREATION ----------------------------------------------------- 	

nameLay = str(tabLayer)
res = nameLay[len(nameLay)-3:len(nameLay)]					## Resolution grid

# Set local variables
cellSize = int(res[0])*100
splineType = "TENSION";
weight = 0.1;

## ----- CHECKERBOARD TEST ------------------------------------ CbT_out_6t1 CbT_out_6t2 CbT_out_18t1 CbT_out_18t1

cbtOutTest = ["cbT_6dT" , "cbT_6daT" , "cbT_15dT" , "cbT_15daT"] 
freqTest = ["6Hz_dT" , "6Hz_daT" , "15Hz_dT" , "15Hz_daT"]


for i in cbtOutTest:
	zFieldIn = i
	outName = "/" + zFieldIn + "_" + res;					## /CbT_out_6t1_(res)km
	outPathName = outFolder + outName;
	outSpline = Spline(nameLay,zFieldIn,cellSize,splineType,weight);
	outSpline.save(outPathName);
	cbtLayer = arcpy.mapping.Layer(outPathName);
	arcpy.mapping.AddLayer(dataFrame,cbtLayer,"BOTTOM");

del outSpline;

for x in cbtOutTest:
	cbtLayer = x + "_" + res
	arcpy.ApplySymbologyFromLayer_management(cbtLayer, "C:/Roby/PhD/5. Deception/GIS/bwColorScale.lyr");



## ----- FREQUENCY MAPS ------------------------------------
	
for j in freqTest:
	zField = j
	outName = "/" + zField + "_" + res;
	outPathName = outFolder + outName;
	outSpline = Spline(nameLay,zField,cellSize,splineType,weight);
	outSpline.save(outPathName);
	newLayer = arcpy.mapping.Layer(outPathName);
	arcpy.mapping.AddLayer(dataFrame,newLayer,"BOTTOM")

del outSpline;
	
for y in freqTest:
	aLayer = y + "_" + res
	arcpy.ApplySymbologyFromLayer_management(aLayer, "C:/Roby/PhD/5. Deception/GIS/magmaColorScale.lyr");
	


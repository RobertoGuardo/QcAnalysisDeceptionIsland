import arcpy, os
from arcpy import env  
from arcpy.sa import *

# Set environment settings
mxd = arcpy.mapping.MapDocument("CURRENT")
dataFrame = arcpy.mapping.MapDocument("CURRENT").activeDataFrame			## Define the dataFrame


## STEP 1: CREATION OF THE TABLE ---------------------------------------- 
table = arcpy.GetParameter(0)									## Load the table
tabName = arcpy.GetParameterAsText(1)							## Set the table name
outFolder = "C:/Roby/PhD/5. Deception/GIS2/" + tabName			## Destination folder Qc_20283_1km
filePath = outFolder + "/" + tabName + ".shp"					## Path to file

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
    arcpy.CreateFolder_management("C:/Roby/PhD/5. Deception/GIS2/", tabName) ## Create destination folder
	
arcpy.FeatureClassToShapefile_conversion(out_Layer, outFolder)

## END STEP 1 ----------------------------------------

tabLayer = arcpy.mapping.Layer(filePath)							## Qc_20283_1km.shp
arcpy.mapping.AddLayer(dataFrame,tabLayer,"BOTTOM")					## Add the layer to the active dataFrame

## STEP 2: MAPS CREATION ----------------------------------------------------- 	

nameLay = str(tabLayer)							## Qc_20283_1km.shp
res = nameLay[len(nameLay)-3:len(nameLay)]		## Resolution grid

sepChar = "_"
iSepChar = nameLay.find(sepChar)


ray = nameLay[iSepChar+1:len(nameLay)-4]

if ray == "20283":
	it = "0"
elif ray == "14972":
	it = "1"
elif ray == "13105":
	it = "2"
elif ray == "7895":
	it = "3"
elif ray == "7196":
	it = "4"

# Set local variables
cellSize = int(res[0])*100
splineType = "TENSION";
weight = 0.1;

## ----- FREQUENCY MAPS

groupLay = "Qc_it" + it + "_" + ray + "_" + res
targetGroupLayer = arcpy.mapping.ListLayers(mxd, groupLay, dataFrame)[0]

for freq in range(0,6):
	zField = str(freq*3+6) + "Hz";  ## To locate the layer in order use str(21-freq*3), otherwise str(freq*3+6)
	outName = "/it"+ it + "_" + zField + "_" + res;
	outPathName = outFolder + outName;
	outSpline = Spline(nameLay, zField, cellSize, splineType, weight);
	outSpline.save(outPathName);
	newLayer = arcpy.mapping.Layer(outPathName);
	arcpy.mapping.AddLayerToGroup(dataFrame, targetGroupLayer, newLayer, "BOTTOM")
	# arcpy.mapping.AddLayer(dataFrame,newLayer,"BOTTOM");
	
for x in range(0,6):
	aLayer = "it" + it + "_" + str(x*3+6) + "Hz" + "_" + res
	arcpy.ApplySymbologyFromLayer_management(aLayer, "C:/Roby/PhD/5. Deception/GIS2/magmaColorScale.lyr");

del outSpline;





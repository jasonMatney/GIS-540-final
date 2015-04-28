'''
Title: Regression Analysis of Spatially Indexed Permafrost Data.
Author: Jason Matney - jamatney
Class: GIS 540 GIS Programming
Created: 4/15/2015
Updated: 4/24/2015
Version: 1.2
Description: A exploratory data analysis and subsequent regression execution is used
to estimate the binomial presence or absence of permafrost within a region of the
Alaskan Discontinuous Zone. This Geo processing solution will incorporate geo referenced
data from the USGS data sets of interest. Four previously derived covariance
will be leveraged for used as predictor variables.
The output in ArcGIS Desktop might be easy for a GIS Analyst to read
and understand, but the non-GIS users, such as drivers, need a simplified
version of the information to do their part. The goal of this tool is to
demonstrate how to generate PDF reports with maps and directions for a Regression
analysis. This solution also takes advantage of ArcGIS Online to store the PDFs
as well as the routes as feature services to use within web maps. The permafrost predictions
can be viewed by analysts in phone applications such as
ArcGIS Online App to view potential research sites and stay on track of Climate Change^TM.
The point is to simplify the process of doing a Regression Analysis and to
produce a user friendly interface for the researcher to follow.
'''


# Import system modules
import arcpy, os
import pandas as pd
from arcpy import env

def find_all_csv(workspace):
    prev_workspace = arcpy.env.workspace

    arcpy.env.workspace = workspace
    csv = arcpy.ListFiles("*csv")

    arcpy.env.workspace = prev_workspace
    return csv


# Begin Script processing
# Permafrost regression processing
# Exploratory Regression of permafrost in the Alaskan Discontinuous Zone
# using the Exploratory Regression Tool


# # User Parameters and variable setup
Workspace = arcpy.GetParameterAsText(0)
coords = pd.read_csv(arcpy.GetParameterAsText(1))
permafrost = pd.read_csv(arcpy.GetParameterAsText(2))
heatload = pd.read_csv(arcpy.GetParameterAsText(3))
temp = pd.read_csv(arcpy.GetParameterAsText(4))
slope = pd.read_csv(arcpy.GetParameterAsText(5))
cti = pd.read_csv(arcpy.GetParameterAsText(6))
texture = pd.read_csv(arcpy.GetParameterAsText(7))
outputfolder = arcpy.GetParameterAsText(8)
shapefile = arcpy.GetParameterAsText(9)


# Set workspace

env.workspace = Workspace

# Set overwrite to TRUE
env.overwriteOutput = True
arcpy.gp.overwriteOutput = True

try:
    env.overwriteOutput = True
    arcpy.AddMessage("Creating spatially referenced csv files...")
    """ cbind coords to covaraites """
    covariates = pd.concat([permafrost, heatload, temp, slope, cti, texture, coords], axis=1)

    # Write to csv
    covariates.to_csv(os.path.join(outputfolder, "covariates.csv"), index=False)
    """Takes a CSV as input, write a SHP as output"""

    arcpy.AddMessage("Converting to shapefile...")
    dir_string = os.path.join(Workspace, outputfolder)
    shp_output_dir = os.path.join(dir_string, "shapefiles")

    files = find_all_csv(os.path.join(Workspace, outputfolder))
    for file in files:
        csv_input = dir_string + r"\{0}".format(file)
        temp_layer = os.path.splitext(os.path.basename(csv_input))[0] # == "input"

        # Get NAD 1983 Alaska Albers Spatial reference
        prjfile = arcpy.SpatialReference(3338)

        arcpy.MakeXYEventLayer_management(csv_input, "long_site", "lat_site", temp_layer, prjfile)
        # exports a feature to a directory as a shapefile.

        # set covaraite shapefile feature class
        fc = os.path.join(shp_output_dir, temp_layer) + ".shp"

        # Delete covariate shapefile if exists
        if arcpy.Exists(fc):
            arcpy.Delete_management(fc)

        arcpy.FeatureClassToShapefile_conversion(temp_layer, shp_output_dir)
        arcpy.Delete_management(temp_layer) # clean up layer, for completeness

    arcpy.AddMessage("Shapefile conversion complete.")


except arcpy.AddError("\tThere may be an issue with your input files."):
    arcpy.AddMessage("\pPlease check covariate csv.")

try:
    arcpy.AddMessage("Starting Permafrost Analysis...")
    # Set geoprocessor object property to overwrite existing output, by default
    arcpy.gp.overwriteOutput = True
    arcpy.env.workspace = os.path.join(Workspace, outputfolder) + "/shapefiles"

    # Define spatial projection
    infc = "covariates.shp"
    # Get NAD 1983 Alaska Albers Spatial reference
    prjfile = arcpy.SpatialReference(3338)
    arcpy.DefineProjection_management(infc, prjfile)

    # create projected outfile
    input_feature = "covariates.shp"
    output_feature_class = "covariates_proj.shp"
    out_coordinate_system = arcpy.SpatialReference(3338)
    arcpy.Project_management(input_feature, output_feature_class, out_coordinate_system)

    # # Join the permafrost feature class to the Covariate feature class
    # # Process: Spatial Join
    fieldMappings = arcpy.FieldMappings()
    fieldMappings.addTable("denali_proj.shp")
    fieldMappings.addTable("covariates_proj_clip.shp")

    arcpy.AddMessage("\tPerforming Spatial Join Analysis on Permafrost Data...")
    sj = arcpy.SpatialJoin_analysis("denali_proj.shp", "covariates_proj_clip.shp", "spatial_join.shp",
                                    "JOIN_ONE_TO_MANY",
                                    "KEEP_ALL",
                                    fieldMappings,
                                    "INTERSECTS", "", "")

    # Delete extra fields to clean up the data
    # Process: Delete Field
    arcpy.AddMessage("\tDeleting Fields...")
    arcpy.DeleteField_management("spatial_join.shp", "SP_ID;OBJECTID;NAME;STATE_NAME;STATE_FIPS;CNTY_FIPS;FIPS;POP2010;POP10_SQMI;POP2013;POP13_SQMI;WHITE;BLACK;AMERI_ES;ASIAN;HAWN_PI;HISPANIC;OTHER;MULT_RACE;MALES;FEMALES;AGE_UNDER5;AGE_5_9;AGE_10_14;AGE_15_19;AGE_20_24;AGE_25_34;AGE_35_44;AGE_45_54;AGE_55_64;AGE_65_74;AGE_75_84;AGE_85_UP;MED_AGE;MED_AGE_M;MED_AGE_F;HOUSEHOLDS;AVE_HH_SZ;HSEHLD_1_M;HSEHLD_1_F;MARHH_CHD;MARHH_NO_C;MHH_CHILD;FHH_CHILD;FAMILIES;AVE_FAM_SZ;HSE_UNITS;VACANT;OWNER_OCC;RENTER_OCC;NO_FARMS12;AVE_SIZE12;CROP_ACR12;AVE_SALE12;SQMI;Shape_Leng;Shape_Area;TARGET_FID;Join_Count;STATE_NAME;DRAWSEQ;STATE_FIPS;SUB_REGION;STATE_ABBR")
    arcpy.AddMessage("\tListing Remaining Fields...")

    layer_fields = arcpy.ListFields("spatial_join.shp")
    for field in layer_fields:
        print "{0} is a type of {1} with a length of {2}"\
            .format(field.name, field.type, field.length)

    # Define spatial projection
    infc = "spatial_join.shp"
    # Get NAD 1983 Alaska Albers Spatial reference
    prjfile = arcpy.SpatialReference(3338)
    arcpy.DefineProjection_management(infc, prjfile)
    sr = arcpy.Describe("spatial_join.shp").spatialReference
    print sr.name

    # create projected outfile
    input_feature = "spatial_join.shp"
    output_feature_class = "spatial_join_proj.shp"
    out_coordinate_system = arcpy.SpatialReference(3338)
    arcpy.Project_management(input_feature, output_feature_class, out_coordinate_system)

    # Create Spatial Weights Matrix for Calculations
    # Process: Generate Spatial Weights Matrix
    arcpy.AddMessage("\tBuilding Spatial Weights Matrix...")
    swm = arcpy.GenerateSpatialWeightsMatrix_stats("spatial_join_proj.shp", "JOIN_FID", "spatial_weights.swm", "K_NEAREST_NEIGHBORS")
    #
    arcpy.AddMessage("\tSpatial Weights Matrix generated...")

    # Exploratory Regression Analysis for permafrost

    arcpy.AddMessage("\tPerforming Exploratory Spatial Regression...")
    er = arcpy.ExploratoryRegression_stats(Input_Features="spatial_join_proj.shp",
                                      Dependent_Variable="permafrost",
                                      Candidate_Explanatory_Variables="heatload;temp;slope;cti;texture",
                                      Weights_Matrix_File="spatial_weights.swm",
                                      Output_Report_File="results.txt",
                                      Maximum_Number_of_Explanatory_Variables="5",
                                      Minimum_Number_of_Explanatory_Variables="1",
                                      Minimum_Acceptable_Adj_R_Squared="0.5",
                                      Maximum_Coefficient_p_value_Cutoff="0.05",
                                      Maximum_VIF_Value_Cutoff="7.5",
                                      Minimum_Acceptable_Jarque_Bera_p_value="0.1",
                                      Minimum_Acceptable_Spatial_Autocorrelation_p_value="0.1")
    arcpy.AddMessage("Regression Analysis Complete.")

except arcpy.ExecuteError("An error occurred during processing:\n"):
    msgs = arcpy.GetMessages(2)
    arcpy.AddError("An error occurred during processing:\n")
    arcpy.AddError(msgs)
    arcpy.AddError("\nP Check that inputs are formatted correctly.")

# # Make RouteDirections object to pull direction information off of it
# d = VRPS.RouteDirection(directions)
#
# # update sublayers name reference
# ordersLayer = arcpy.mapping.ListLayers(mxd, "Orders")[0]
# depotsLayer = arcpy.mapping.ListLayers(mxd, "Depots")[0]
# routesLayer = arcpy.mapping.ListLayers(mxd, "Routes")[0]
#
#
# # Start mapbook and upload processing for each inspector
# arcpy.AddMessage("Starting Mapbook processing...")
# routesCursor = arcpy.da.SearchCursor(routestable, ["Name"])
# routebookcollection = []
# for inspector_row in routesCursor:
#     Name = inspector_row[0]
#     # Create empty Route book PDF
#     arcpy.AddMessage("\tStarting {}'s mapbook...".format(Name))
#     pdf_path = os.path.join(outputfolder, "{0}_RouteBook_{1}.pdf".format(Name, date))
#     pdf = arcpy.mapping.PDFDocumentCreate(pdf_path)
#     routebookcollection.append(pdf_path)
#     # select individual inspector orders
#     arcpy.SelectLayerByAttribute_management(depotsLayer, "NEW_SELECTION",\
#             "Name = 'Assessors Office'")
#     arcpy.SelectLayerByAttribute_management(ordersLayer, "NEW_SELECTION",\
#             "RouteName = '{}'".format(Name))
#     count = int(arcpy.GetCount_management(ordersLayer).getOutput(0))
#     sequence_num = 2
#     # Create a temporary folder in the output to build order pages
#     outputfolder_temp = os.path.join(outputfolder, Name + "_temp")
#     if arcpy.Exists(outputfolder_temp):
#         arcpy.Delete_management(outputfolder_temp)
#     os.makedirs(outputfolder_temp)
#     # start page build
#     while sequence_num <= (count + 1):
#         if sequence_num == 2:
#             arcpy.SelectLayerByAttribute_management(ordersLayer, \
#                     "NEW_SELECTION", \
#                     "Sequence = 2 AND RouteName = '{0}'".format(Name))
#         else:
#             arcpy.SelectLayerByAttribute_management(depotsLayer, \
#                         "CLEAR_SELECTION")
#             arcpy.SelectLayerByAttribute_management(ordersLayer, \
#                 "NEW_SELECTION", \
#                 "(Sequence = {0} OR Sequence = {1}) AND RouteName = '{2}'".format(\
#                 sequence_num, (sequence_num - 1), Name))
#         df.zoomToSelectedFeatures()
#         df.scale = df.scale + 2000
#         # Find directions for inspector for select order and update map
#         d.findStringPositions(Name)
#         dic = d.seekLines()
#         DirecttextElement = arcpy.mapping.ListLayoutElements(mxd, \
#                          "TEXT_ELEMENT", "directions")[0]
#         DirecttextElement.text = dic[(sequence_num - 1)]
#         DirecttextElement.elementWidth = 3.25
#         InspectTexttElement = arcpy.mapping.ListLayoutElements(mxd, \
#                          "TEXT_ELEMENT", "inspection")[0]
#         ordersCursor = arcpy.da.SearchCursor(ordersLayer, ["Name"])
#         for row in ordersCursor:
#             InspectTexttElement.text = row[0]
#         page_name = os.path.join(outputfolder_temp, "{0}_{1}.pdf".format(Name,\
#                                  sequence_num))
#         # Export map and apped it to main route book pdf
#         arcpy.mapping.ExportToPDF(mxd, page_name, "PAGE_LAYOUT")
#         pdf.appendPages(page_name)
#         sequence_num += 1
#     pdf.saveAndClose()
#     arcpy.Delete_management(outputfolder_temp)
#     arcpy.AddMessage("\t{}'s Mapbook created.".format(Name))
#
#     # upload routes to Agol and publish
#     arcpy.AddMessage("\tStarting upload of Route and Orders shapefiles...")
#     upload_routes = VRPS.uploadPublish(Name, date, outputfolder, \
#                     routesLayer, "Name = '{}'".format(Name), username, token)
#     upload_orders = VRPS.uploadPublish(Name, date, outputfolder, \
#                     ordersLayer, "RouteName = '{}'".format(Name), username, token)
#     if upload_orders != None and upload_routes != []:
#         VRPS.makeWebmap(Name, date, upload_routes, upload_orders, username, token)
#     arcpy.AddMessage("\tFinsihed processing {}'s pdf and webmap.\n".format(Name))
#
#
# # Upload PDF to ArcGIS Online and share them with the organization
# arcpy.AddMessage("Starting PDF upload process...")
# output_pdf_dict = {}
# for routebook in routebookcollection:
#     upload_pdf_dict = VRPS.uploadPDF(routebook, username, token)
#     output_pdf_dict.update(upload_pdf_dict)
# VRPS.sharePDFs(output_pdf_dict, username, token)
#
# # final cleanup
# del d, inspector_row, routesCursor, mxd, df, vprLayer, ordersLayer
# del routesLayer, depotsLayer

# arcpy.AddMessage("End of Processing")

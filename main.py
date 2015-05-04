"""
Title: Regression Analysis of Spatially Indexed Permafrost Data.
Author: Jason Matney - jamatney
Class: GIS 540 GIS Programming
Created: 4/15/2015
Updated: 4/30/2015

Description: A exploratory data analysis and subsequent regression execution is used
to estimate the prowess of predictive variables behind binomial presence or absence of permafrost
within a region of the Alaskan Discontinuous Zone.
This Geo processing solution will incorporate geo referenced
data from the USGS data sets of interest. Four previously derived covariance
will be leveraged for used as predictor variables.

The goal of this tool is to demonstrate how to generate mxd reports
with maps for a Regression analysis.

Instructions
1) use attached rpy2Pandas.pdf to install rpy and pandas
2) When R prompts you to use a personal directory, select Yes,
then select a CRAN mirror nearest to you (USA-MD works well)
3) If Exploratory Regression Analysis fails, Comment out 289 - 308 and
uncomment 278 - 284 to use GWR for the purpose of having output.
"""

# Import system modules
import arcpy, os
import pandas as pd
from arcpy import env
import rpy2.robjects as ro
import shutil

class SpatialJoin:
    def __init__(self, feature_shp, covariate_shp, outfile_shp, spatial_ref):
        self.feature_shp = feature_shp
        self.covariate_shp = covariate_shp
        self.outfile_shp = outfile_shp
        self.spatial_ref = spatial_ref
        self.output_feature_class = None

    def spatial_join(self):
        fieldMappings = arcpy.FieldMappings()
        fieldMappings.addTable(self.feature_shp)
        fieldMappings.addTable(self.covariate_shp)

        arcpy.AddMessage("\tPerforming Spatial Join Analysis on Permafrost Data...")
        arcpy.SpatialJoin_analysis(target_features=self.feature_shp,
                                   join_features=self.covariate_shp,
                                   out_feature_class=self.outfile_shp,
                                   join_operation="JOIN_ONE_TO_MANY",
                                   field_mapping=fieldMappings,
                                   match_option="COMPLETELY_CONTAINS")

        arcpy.AddMessage("\tOutput is 'spatial_join.shp'...")
        # Delete extra fields to clean up the data
        # Process: Delete Field
        arcpy.AddMessage("\tDeleting Fields...")
        arcpy.DeleteField_management(self.outfile_shp,   "NAME;STATE_NAME;"
                                                         "STATE_FIPS;CNTY_FIPS;"
                                                         "FIPS;POP2010;POP10_SQMI;"
                                                         "POP2013;POP13_SQMI;WHITE;"
                                                         "BLACK;AMERI_ES;ASIAN;HAWN_PI;"
                                                         "HISPANIC;OTHER;MULT_RACE;MALES;"
                                                         "FEMALES;AGE_UNDER5;AGE_5_9;"
                                                         "AGE_10_14;AGE_15_19;AGE_20_24;"
                                                         "AGE_25_34;AGE_35_44;AGE_45_54;"
                                                         "AGE_55_64;AGE_65_74;AGE_75_84;"
                                                         "AGE_85_UP;MED_AGE;MED_AGE_M;"
                                                         "MED_AGE_F;HOUSEHOLDS;AVE_HH_SZ;"
                                                         "HSEHLD_1_M;HSEHLD_1_F;MARHH_CHD;"
                                                         "MARHH_NO_C;MHH_CHILD;FHH_CHILD;"
                                                         "FAMILIES;AVE_FAM_SZ;HSE_UNITS;"
                                                         "VACANT;OWNER_OCC;RENTER_OCC;NO_FARMS12;"
                                                         "AVE_SIZE12;CROP_ACR12;AVE_SALE12;"
                                                         "SQMI;Shape_Leng;Shape_Area;"
                                                         "Join_Count;STATE_NAME;DRAWSEQ;STATE_FIPS;"
                                                         "SUB_REGION;STATE_ABBR")
        arcpy.AddMessage("\tListing Remaining Fields...")

        layer_fields = arcpy.ListFields(self.outfile_shp)
        for field in layer_fields:
            print "{0} is a type of {1} with a length of {2}"\
                .format(field.name, field.type, field.length)

    def define_and_project(self, ofc):
        self.output_feature_class = ofc
        """Define spatial projection and create projected outfile"""
        prjfile = arcpy.SpatialReference(self.spatial_ref)
        arcpy.DefineProjection_management(self.outfile_shp, prjfile)
        sr = arcpy.Describe(self.outfile_shp).spatialReference
        print "Your ouput feature calss is : {0}".format(self.output_feature_class)
        print "Your spatial reference code corresponds to: {0}".format(sr.name)

        # Create projected outfile
        arcpy.Project_management(self.outfile_shp, self.output_feature_class, self.spatial_ref)


def find_all_csv(workspace):
    prev_workspace = arcpy.env.workspace

    arcpy.env.workspace = workspace
    csv = arcpy.ListFiles("*csv")

    arcpy.env.workspace = prev_workspace
    return csv


def shapefile_conversion(file_list, csv_location, output_location, spatial_reference_code):
        """Takes a CSV as input, write a SHP as output"""
        arcpy.AddMessage("Converting to shapefile...")
        for file in file_list:
            csv_input = csv_location + r"\{0}".format(file)
            temp_layer = os.path.splitext(os.path.basename(csv_input))[0] # == "input"
            # Get NAD 1983 Alaska Albers Spatial reference
            prjfile = arcpy.SpatialReference(spatial_reference_code)
            arcpy.MakeXYEventLayer_management(csv_input, "long_site", "lat_site", temp_layer, prjfile)
            # exports a feature to a directory as a shapefile.
            # set covaraite shapefile feature class
            fc = os.path.join(output_location, temp_layer) + ".shp"
            # Delete covariate shapefile if exists
            if arcpy.Exists(fc):
                arcpy.Delete_management(fc)
            arcpy.FeatureClassToShapefile_conversion(temp_layer, output_location)
            arcpy.Delete_management(temp_layer) # clean up layer, for completeness
        arcpy.AddMessage("Shapefile conversion complete.")


def clip_features(inf, clip, out):
    arcpy.AddMessage("\tClipping data to Denali Counties...")

    # Execute Clip
    arcpy.Clip_analysis(inf, clip, out)
    arcpy.AddMessage("\tClip Executed, Leveraging for further analysis...")


def make_map(dir):
    prev_workspace = arcpy.env.workspace
    arcpy.env.workspace = dir
    #Make a Map Document
    try:
        mapList = arcpy.ListFiles("*.mxd")
        mxd_path = os.path.join(dir, mapList[0])

        mxd = arcpy.mapping.MapDocument(mxd_path)
        data_frames = arcpy.mapping.ListDataFrames(mxd)
        data_frame = data_frames[0]
        fcs = arcpy.ListFeatureClasses()

        for f in fcs:
            out_layer = f[:-4]
            out_layer_file = out_layer + ".lyr"
            arcpy.MakeFeatureLayer_management(f, out_layer)
            arcpy.SaveToLayerFile_management(out_layer, out_layer_file)
            layer_object = arcpy.mapping.Layer(f)
            arcpy.mapping.AddLayer(data_frame, layer_object)

        arcpy.RefreshTOC()
        project_map = map_document[:-4] + "Presentation.mxd"
        mxd.saveACopy(os.path.join(prev_workspace, project_map))
        os.startfile(os.path.join(prev_workspace, project_map))
        arcpy.env.workspace = prev_workspace
        del mxd
    except arcpy.AddMessage("\tPlease check inputs"):
        arcpy.GetMessages()



# User Parameters and variable setup
Workspace = arcpy.GetParameterAsText(0)  # C:\Users\jamatney\Desktop\GIS540_project_jamatney

# load covariates
coords = pd.read_csv(arcpy.GetParameterAsText(1))  # data\MegaAlaska\covariates\coords.csv
permafrost = pd.read_csv(arcpy.GetParameterAsText(2))  # data\MegaAlaska\covariates\permafrost.csv
heatload = pd.read_csv(arcpy.GetParameterAsText(3))  # data\MegaAlaska\covariates\heatload.csv
temp = pd.read_csv(arcpy.GetParameterAsText(4))  # data\MegaAlaska\covariates\temp.csv
slope = pd.read_csv(arcpy.GetParameterAsText(5))  # data\MegaAlaska\covariates\slope.csv
cti = pd.read_csv(arcpy.GetParameterAsText(6))  # data\MegaAlaska\covariates\cti.csv
texture = pd.read_csv(arcpy.GetParameterAsText(7))  # data\MegaAlaska\covariates\texture.csv

# set output folders
outputfolder = arcpy.GetParameterAsText(8)  # data
denali = arcpy.GetParameterAsText(9)  # input\discontinuous_zone.shp
projcode = int(arcpy.GetParameterAsText(10))  # 3338
projected_covariates = arcpy.GetParameterAsText(11)  # covariates_proj.shp
join_file = arcpy.GetParameterAsText(12)  # spatial_join.shp
map_document = arcpy.GetParameterAsText(13)  # project_jamatney.mxd
denali_cov = os.path.basename(denali)[:-4] + "_covariates.shp"

# Environmental Variables
arcpy.env.workspace = Workspace
arcpy.env.overwriteOutput = True
arcpy.gp.overwriteOutput = True

# Generated initial variables
dir_string = os.path.join(Workspace, outputfolder)
shp_output_dir = os.path.join(dir_string, "shapefiles")
files = find_all_csv(dir_string)
print arcpy.env.workspace
csv_files = find_all_csv(os.path.join(dir_string, "MegaAlaska\covariates"))
print csv_files
for f in csv_files:
    print "Leveraging Covariate file: {0}".format(f)

#############
# Data Prep #
#############

try:
    env.overwriteOutput = True
    arcpy.AddMessage("Creating spatially referenced csv files...")
    """ cbind coords to covaraites """

    covariates = pd.concat([permafrost, heatload, temp, slope, cti, texture, coords], axis=1)

    # Write to csv
    covariates.to_csv(os.path.join(outputfolder, "covariates.csv"), index=False)
    # Convert ot Shapefile
    shapefile_conversion(files, dir_string, shp_output_dir, projcode)

except arcpy.AddError("\tThere may be an issue with your input files."):
    arcpy.AddMessage("\pPlease check covariate csv.")

try:
    arcpy.AddMessage("Starting Permafrost Analysis...")
    # Set geoprocessor object property to overwrite existing output, by default
    arcpy.gp.overwriteOutput = True
    arcpy.env.workspace = shp_output_dir

    arcpy.AddMessage("\tReprojecting data using R...")

    # Reprojection done in R
    ro.globalenv['dsn'] = arcpy.env.workspace
    ro.r('''
    if (!require("pacman")) install.packages("pacman")
    pacman::p_load(sp, rgdal, raster)
    ak_proj <- readOGR(dsn, "Alaska")
    covariates  <- readOGR(dsn,"covariates")
    cov.a  <- project(covariates@coords, proj4string(ak_proj))
    cov.b <- cbind(cov.a, covariates@data[,1:6])
    colnames(cov.b) <- c("x","y","permafrost","heatload","temp","slope","cti","texture")
    coordinates(cov.b) <-~x+y
    proj4string(cov.b) <- proj4string(ak_proj)
    cov.ab <- spTransform(cov.b,CRS(proj4string(ak_proj)))
    ''')
    arcpy.Delete_management("covariates.shp")
    arcpy.Delete_management("covariates.lyr")
    if arcpy.Exists(projected_covariates):
            arcpy.Delete_management(projected_covariates)
    ro.r('writeOGR(cov.ab, dsn, layer="covariates_proj", driver="ESRI Shapefile")')

except arcpy.AddMessage("An error occurred during processing:\n"):
    arcpy.GetMessages()

####################
# Spatial Analysis #
####################

# Change working directory
arcpy.env.workspace = shp_output_dir
arcpy.gp.overwriteOutput = True

# Hardcode generic spatial join file names
denali_base = os.path.basename(denali)

# Perform Clip
clip_features(projected_covariates, denali, denali_cov)

# Join the permafrost feature class to the covariate feature class
# Process: Spatial Join
myJoin = SpatialJoin(denali, denali_cov, join_file, projcode)
myJoin.spatial_join()

try:
    # If exploratory regression fails use GWR
    # arcpy.AddMessage("\tGeographically Weighted Regression...")
    # Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
    # The following inputs are layers or table views: "spatial_join"
    # gwr = arcpy.GeographicallyWeightedRegression_stats(in_features=join_file,dependent_field="permafrost",
    #                                                    explanatory_field="heatload;temp;slope;cti",
    #                                                    out_featureclass="GeographicallyWeightedRegression",
    #                                                    kernel_type="FIXED",bandwidth_method="AICc",
    #                                                    distance="#",number_of_neighbors="30",
    #                                                    weight_field="#",coefficient_raster_workspace="#",
    #                                                    cell_size="0.306677063650721",in_prediction_locations="#",
    #                                                    prediction_explanatory_field="#",out_prediction_featureclass="#")

    #Create Spatial Weights Matrix for Calculations
    #Process: Generate Spatial Weights Matrix
    arcpy.AddMessage("\tBuilding Spatial Weights Matrix...")
    swm = arcpy.GenerateSpatialWeightsMatrix_stats(Input_Feature_Class=join_file,
                                                   Unique_ID_Field="JOIN_FID",
                                                   Output_Spatial_Weights_Matrix_File="spatial_weights.swm",
                                                   Conceptualization_of_Spatial_Relationships="K_NEAREST_NEIGHBORS")
    arcpy.AddMessage("\tSpatial Weights Matrix generated...")

    #Exploratory Regression Analysis for permafrost
    arcpy.AddMessage("\tPerforming Exploratory Spatial Regression...")
    er = arcpy.ExploratoryRegression_stats(Input_Features=join_file,
                                           Dependent_Variable="permafrost",
                                           Candidate_Explanatory_Variables=\
                                           "heatload;temp;slope;cti;texture",
                                           Weights_Matrix_File="spatial_weights.swm",
                                           Output_Report_File="results.txt",
                                           Maximum_Number_of_Explanatory_Variables="5",
                                           Minimum_Number_of_Explanatory_Variables="1",
                                           Minimum_Acceptable_Adj_R_Squared="0.3",
                                           Maximum_Coefficient_p_value_Cutoff="0.10",
                                           Maximum_VIF_Value_Cutoff="7.5")

    # Move results.txt to data folder
    if arcpy.Exists( os.path.join(shp_output_dir, "results.txt")):
        shutil.move( os.path.join(shp_output_dir, "results.txt"),
                     os.path.join(dir_string, "regression_results\\results.txt"))

    if arcpy.Exists( os.path.join(shp_output_dir, "results.txt.xml")):
        shutil.move( os.path.join(shp_output_dir, "results.txt.xml"),
                     os.path.join(dir_string, "regression_results\\results.txt.xml"))


    if arcpy.Exists(join_file):
        arcpy.Delete_management(join_file)

    arcpy.AddMessage("Geographically Weighted Regression Analysis Complete.")

except arcpy.ExecuteError("An error occurred during processing:\n"):
    msgs = arcpy.GetMessages(2)
    arcpy.AddError(msgs)
    arcpy.AddError("\nP Check that inputs are formatted correctly.")

# Reset workspace to home directory
arcpy.env.workspace = Workspace

# Make a Map
make_map(shp_output_dir)

arcpy.AddMessage("End of Processing.")
arcpy.AddMessage("Please find Exploratory Regression Results in data/regression_results/results.txt")

arcpy.AddMessage("******* Output from R Console Below *******")

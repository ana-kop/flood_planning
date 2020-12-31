import rasterio.plot
import rasterio.mask
from rasterio.windows import Window
from rasterio.plot import show
from rasterio.plot import plotting_extent
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, Polygon, box, LineString
from shapely.ops import nearest_points
import pyproj
import numpy as np

osgb36 = pyproj.Proj('EPSG:27700')      # Set CRS to BNG

# Importing data
path = "/Users/danniharnett/Desktop/Material"
background = rasterio.open(path+"/background/raster-50k_2724246.tif")
elevation = rasterio.open(path+"/elevation/SZ.asc", "r")


# # ITN data
# itn = pd.read_json(path + "itn/solent_itn.json")
# # print(itn.head())
# # print(list(itn))
# #
# # # Nodes and links data
# nodes = gpd.read_file(path + "/roads/nodes.shp")
# # # print(nodes.head())
# links = gpd.read_file(path + "/roads/links.shp")
# # # print(links.head())


'''  Task 1:
Ask the user to input their current location as a British National Grid coordinate (easting and northing),
test whether the user is within a box (430000, 80000) and (465000, 95000); quit program if not
'''
#
print(" \n This program will help you find the the quickest route to walk to the highest point of land \
within a 5km radius.")
#
try:
   x = float(input("Please enter the Easting coordinate of your location: "))
# If user entered not a number and there is ValueError, then:
except ValueError:
   x = float(input("Please input a NUMBER. Please enter the X-coordinate of your point: "))

try:
   y = float(input("Please enter the Northing coordinate of your location: "))
# If user entered not a number and there is ValueError, then:
except ValueError:
   y = float(input("Please input a NUMBER. Please enter the Y-coordinate of your point: "))

# test whether the user is within a box (430000, 80000) and (465000, 95000).
inside = True

# code for testing that user is within the (430000, 80000) and (465000, 95000)

if inside is False:
   print("Unable to assist in finding highest point of land.")


location = Point(x, y)  # User location
buf = location.buffer(5000)             # Create 5km buffer polygon.


''' TASK 2: 
 Find highest point within 5km radius
'''

buffer_gdf = gpd.GeoDataFrame({'geometry': buf}, index=[0], crs=osgb36)  # Put buffer polygon into geodataframe (gdf)

# Mask elevation with buffer polygon & crop. Points outside of the buffer set to value of -100
elevation_output_image, output_transform = rasterio.mask.mask(elevation, buffer_gdf['geometry'],
                                                              nodata=-100, crop=True, filled=True)
out_meta = elevation.meta       # match the output image metadata with elevation metadata
out_meta.update({"driver": "GTiff",     # update metadata for elevation output image
                 "height": elevation_output_image.shape[1],
                 "width": elevation_output_image.shape[2],
                 "transform": output_transform})

with rasterio.open(path+"/elevation_output.tif", "w", **out_meta) as dest:   # write file
    dest.write(elevation_output_image)

elevation.close()

radius = rasterio.open(path+"/elevation_output.tif", "r")       # 5km radius image file opened in read mode
radius_array = radius.read(1)                                    # radius read as numpy array

max_height = np.max(radius_array)   # Find max value within buffer
pix_location_y, pix_location_x = np.where(radius_array == max_height)   # Finds row/columns of pixel value
highest_points = []


def test_transformation(point):  # Ensure transformed coordinates fall within the bounds of the elevation file (BNG)
    if point.x < 42500 or point.x > 470000 or point.y < 75000 or point.y > 100000:
        return True


for i in range(len(pix_location_x)):
    bng_pixel_location = radius.transform*(pix_location_x[i], pix_location_y[i])  # Transforms (row,column) into (x,y)
    high_point = Point(bng_pixel_location)  # Create shapely point at highest point
    if test_transformation(high_point):
        raise Exception("Transformation Error: The determined coordinates do not fall within the bounds of the input \
         file. Please ensure all files use the British National Grid CRS. ")
    else:
        highest_points.append(high_point)



#
# ''' Task 3:
# Identify the nearest ITN node to the user and the nearest ITN node to the highest point identified in Task 2.
# '''

# For every line in nodes, take point from geometry, append point to list, pass it to Multipoint
# https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
# nodes_list = []
# for i in range(len(nodes)):
#     point = nodes.at[i, 'geometry']
#     nodes_list.append(point)
# node_locations = MultiPoint(nodes_list)
#
# nearest_to_user = nearest_points(location, node_locations)
# # print(nearest_geoms[0]) - original point
# nearest_node_user = nearest_to_user[1]
# #print(nearest_node_user)
#
# nearest_to_highground = nearest_points(highest_points[0], node_locations)
# nearest_node_highground = nearest_to_highground[1]
# #print(nearest_node_highground)
#
#
# print('Please proceed to highground at '
#       + str(nearest_node_highground) + '. The closest road for you to go there starts at point '
#       + str(nearest_node_user) + '.')



#  TASK 4:





''' Task 5: 
Plot the user’s starting point, the highest point within the buffer and the shortest route
calculated. Add a color-bar, north arrow, scale bar and a legend
'''

minx = location.x - 5000  # Bounding box dimensions for plotting the 10kmx10km background file
maxx = location.x + 5000
miny = location.y - 5000
maxy = location.y + 5000
bbox = box(minx, miny, maxx, maxy)

bbox_gdf = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs=osgb36)  # Puts bounding box into geodataframe

# crop the background image to 10km x 10km using bounding box:
cropped_background_img, background_transform = rasterio.mask.mask(background, bbox_gdf['geometry'], crop=True)
background_meta = background.meta       # match cropped background img metadata with original background metadata

background_meta.update({"driver": "GTiff",           # update metadata for output img
                        "height": cropped_background_img.shape[1],
                        "width": cropped_background_img.shape[2],
                        "transform": background_transform})

with rasterio.open(path+"/background_output.tif", "w", **background_meta) as dest:  # write file
    dest.write(cropped_background_img)

background.close()

cropped_background = rasterio.open(path+"/background_output.tif", "r")   # opens cropped background
background_extent = plotting_extent(cropped_background)  # background_extent used for plotting


circle = np.ma.masked_where(radius_array == -100, radius_array)  # MASKING ELEVATION FILE TO PLOT CIRCULAR BUFFER

# Masking background sea level values to re-colour
removed_sea = np.ma.masked_where(cropped_background.read(1) == 65, cropped_background.read(1))

# PLOT EVERYTHING:
fig, ax = plt.subplots(figsize=(8, 6))
ax.set_facecolor("paleturquoise")       # Sets sea colour
ax.imshow(removed_sea, cmap="terrain", extent=background_extent)  # Plots cropped background
elevation_radius = ax.imshow(circle, alpha=0.55, aspect=1, extent=background_extent, cmap="viridis")  # Plots elevation
ax.plot(location.x, location.y, 'o', color='tomato', markersize=8, label="User Location")
ax.plot(highest_points[0].x, highest_points[0].y, '^', color='black', markersize=8, label="Highest Point")
cbar = plt.colorbar(elevation_radius, orientation='vertical', fraction=0.025, pad=0.12)
cbar.set_label('Elevation (m)')
plt.arrow((minx+500), (maxy-1200), 0, 950, width=120, head_length=300, length_includes_head=True, facecolor="black",
          edgecolor="black")   # Plot north arrow
ax.plot([(maxx-2500), (maxx-500)], [(miny+475), (miny+475)], color="black", linewidth=1.5)  # Plot 2km scale bar:
plt.text(maxx-1750, miny+125, "2km", fontsize=8, fontweight='bold')   # Plot scale bar label
plt.legend(bbox_to_anchor=(1.025, 1), loc='upper left')
plt.tight_layout()  # Fit everything inside the figure
plt.rcParams.update({'figure.autolayout': True})
plt.show()

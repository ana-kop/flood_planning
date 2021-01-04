import json
import networkx as nx
import rasterio.mask
from rasterio.plot import plotting_extent
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, Polygon, box, LineString, MultiLineString
from shapely.ops import linemerge
import shapely.ops
import numpy as np
from pyproj import Transformer

"""Import data"""
# Elevation data
path = "/Users/ak/Desktop/a_2/Material/"
background = rasterio.open(path + "/background/raster-50k_2724246.tif")
elevation = rasterio.open(path + "/elevation/SZ.asc")
# count = 1 // CRS = BNG
# rasterio.plot.show(elevation)

# Nodes and links data
nodes = gpd.read_file(path + "roads/nodes.shp")
# print(nodes.head())
# print(list(nodes))
# with pandas.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified
#      print(nodes)

links = gpd.read_file(path + "roads/links.shp")

"""  Task 1:
Ask the user to input their current location as a British National Grid coordinate (easting and northing),
test whether the user is within a box (430000, 80000) and (465000, 9500); quit program if not
"""

# print("This program will help you find the the quickest route to walk to the highest point of land "
#      "within a 5km radius from your current location.")
#
# try:
#    x = float(input("Please enter the Easting coordinate of your location: "))
# # If user entered not a number and there is ValueError, then:
# except ValueError:
#    x = float(input("Please input a NUMBER. Please enter the X-coordinate of your point: "))
#
# try:
#    y = float(input("Please enter the Northing coordinate of your location: "))
# # If user entered not a number and there is ValueError, then:
# except ValueError:
#    y = float(input("Please input a NUMBER. Please enter the Y-coordinate of your point: "))
#
# test whether the user is within a box (430000, 80000) and (465000, 95000) - from Yolanda
# if x <= 430000 or x >= 465000 or y <= 80000 or y >= 95000:
#     print("Unable to assist in finding highest point of land.")
#     exit()
# else:
#     print("Thanks for inputting")

# x = 449923.625
# y = 89243.008

x = 439619
y = 85800

# x = 433600 # - problem
# y = 86600

# x = 449000
# y = 96000

#bottom = 451700, 76700
#top = 449000, 96000
#right = 465500, 87700
#left = 433600, 86600
#aldo's = 439619, 85800

user_location = Point(x, y)  # Example user location
buf = user_location.buffer(5000)  # create 5km buffer polygon.


"""Task 2: Find highest point within 5km radius
"""
# https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html

buffer_gdf = gpd.GeoDataFrame({'geometry': buf}, index=[0], crs="EPSG:27700")  # Put buffer polygon into geodataframe

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
        raise Exception("Transformation Error: The determined coordinates do not fall within the bounds of the input "
                        "file. Please ensure all files use the British National Grid CRS. ")
    else:
        highest_points.append(high_point)

"""Task 3.
Identify the nearest ITN node to the user and the nearest ITN node to the highest point identified in Task 2.
"""

# https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
# For every node in the nodes gdf, take the corresponding point from the geometry column and append it to nodes_list.
nodes_list = []
for i in range(len(nodes)):
    node_point = nodes.at[i, 'geometry']
    if node_point.distance(user_location) < 5000:
        nodes_list.append(node_point)
    else:
        pass
nodes_locations = MultiPoint(nodes_list)  # pass the list of node points to a Multipoint shapely object

# node_point = nodes.at[5, 'geometry']
# print(node_point.distance(user_location))

def nearest_itn(point_loc, nodes_loc):
    """Definition of the function that returns the ITN closest to a given point."""
    # Given a location of interest and a Multipoint of ITN nodes, find the nearest ITN node using the nearest_points()
    # method from shapely.ops.
    nearest_nodes_tup = shapely.ops.nearest_points(point_loc, nodes_loc)
    # Select the first ITN node from the resulting tuple.
    nearest_node = nearest_nodes_tup[1]  # nearest_nodes_tup[0] - original point
    return nearest_node


nearest_node_user = nearest_itn(user_location, nodes_locations)
nearest_node_highground = nearest_itn(highest_points[0], nodes_locations)

# Find the FID of a given ITN node
# https://www.geeksforgeeks.org/different-ways-to-iterate-over-rows-in-pandas-dataframe/

# Make a dictionary that allows mapping between the coordinate location of an ITN and its FID
itn_fid_dict = {}
for i in nodes.index:
    geom = nodes['geometry'][i]
    fid_ = nodes['fid'][i]
    geom = str(geom.x) + ',' + str(geom.y)
    itn_fid_dict[geom] = fid_


def coord_to_fid(itn_coord):
    """Definition of the function that takes coordinates of a given ITN node and returns its FID"""
    point_to_convert = str(itn_coord.x) + ',' + str(itn_coord.y)
    fid = itn_fid_dict[point_to_convert]
    return fid


nearest_node_user = coord_to_fid(nearest_node_user)  # POINT (449923.625 89243.008) osgb4000000026146674
nearest_node_highground = coord_to_fid(nearest_node_highground)  # POINT (448627 88053) osgb4000000026145943

""" Task 4: Shortest path.
json structure
{
"roadlinks":
    {"osgb4000000026157611":
        {"length": 358.58, "coords": [   ], "start": "osgb4000000026219225","end": "osgb4000000026141678",
            "natureOfRoad": "Single Carriageway", "descriptiveTerm": "Private Road - Restricted Access"},
}
"""

# Create an empty MultiDiGraph structure
G = nx.MultiDiGraph()

# Load the ITN data
with open(path + "itn/solent_itn.json") as f:
    road_links = json.load(f)['roadlinks']

# Adding edges with weights to the graph G
for i, link_fid in enumerate(road_links):
    length = road_links[link_fid]['length']
    road_coords = road_links[link_fid]['coords']

    road_coords_elev = []
    for i, coord in enumerate(road_coords):
        x = coord[0]
        y = coord[1]
        elev = list(elevation.sample([(x, y)]))[0][0]
        road_coords_elev.append(elev)

    elev = 0
    for i in range(len(road_coords_elev)-1):
        j = i + 1
        b = road_coords_elev[j]
        a = road_coords_elev[i]
        if b > a:
            elev = elev + (b - a)
        else:
            pass

    # Formula to calculate the weight of the node - i.e. the time required to cover it according to the Naismith's rule.
    time_to_walk = ((3 * length) / 250) + (elev / 10)

    # Adding edges with their corresponding weights
    G.add_edge(road_links[link_fid]['start'], road_links[link_fid]['end'], fid=link_fid, weight=time_to_walk)

    # Calculate weights and add edges in the opposite direction
    road_coords_elev.reverse()
    for i in range(len(road_coords_elev)-1):
        j = i + 1
        b = road_coords_elev[j]
        a = road_coords_elev[i]
        if b > a:
            elev = elev + (b - a)
        else:
            pass

    time_to_walk = ((3 * length) / 250) + (elev/10)
    G.add_edge(road_links[link_fid]['end'], road_links[link_fid]['start'], fid=link_fid, weight=time_to_walk)

# Finding the shortest path using the dijkstra algorithm from networkx.
try:
    shortest_path = nx.dijkstra_path(G, nearest_node_user, nearest_node_highground, weight="weight")
except nx.NetworkXNoPath:
    print('Unable to assist in finding the shortest path.')
    exit()

# Converting the FIDs returned by the algorithm into a linestring for plotting; getting path length.
# https://gis.stackexchange.com/questions/223447/weld-individual-line-segments-into-one-linestring-using-shapely
shortest_path_lines = []
path_length = 0
for i, c in enumerate(shortest_path[:-1]):
    fid = G.get_edge_data(c, shortest_path[i + 1])[0]['fid']
    coords = road_links[fid]['coords']
    shortest_path_lines.append(LineString(coords))
    path_length = path_length + road_links[fid]['length']
shortest_path_final = linemerge(MultiLineString(shortest_path_lines))
path_length = round(path_length/1000, 3)

"""Task 5: Output
"""
# bounding box dimensions (for plotting the background file):
minx = user_location.x - 5000  # 5km in each direction
maxx = user_location.x + 5000
miny = user_location.y - 5000
maxy = user_location.y + 5000
bbox = box(minx, miny, maxx, maxy)

# Bounding box dimensions for plotting the 10kmx10km background file
bbox_gdf = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs="EPSG:27700")

# Crop the background image to 10km x 10km using bounding box:
cropped_background_img, background_transform = rasterio.mask.mask(background, bbox_gdf['geometry'], crop=True)
background_meta = background.meta       # Match cropped background img metadata with original background metadata

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
ax.set_facecolor("paleturquoise")  # Sets sea colour
ax.imshow(removed_sea, cmap="terrain", extent=background_extent)  # Plots cropped background
elevation_radius = ax.imshow(circle, alpha=0.55, aspect=1, extent=background_extent, cmap="viridis")  # Plots elevation
ax.plot(user_location.x, user_location.y, 'o', color='tomato', markersize=8, label="User Location")
ax.plot(highest_points[0].x, highest_points[0].y, '^', color='black', markersize=8, label="Highest Point")
cbar = plt.colorbar(elevation_radius, orientation='vertical', fraction=0.025, pad=0.12)
cbar.set_label('Elevation (m)')
plt.arrow((minx+500), (maxy-1200), 0, 950, width=120, head_length=300, length_includes_head=True, facecolor="black",
          edgecolor="black")   # Plot north arrow
ax.plot([(maxx-2500), (maxx-500)], [(miny+475), (miny+475)], color="black", linewidth=1.5)  # Plot 2km scale bar:
plt.text(maxx-1750, miny+125, "2km", fontsize=8, fontweight='bold')   # Plot scale bar label

x, y = shortest_path_final.xy
ax.plot(x, y, color='black', linewidth=2, zorder=1, label="shortest path")

plt.legend(bbox_to_anchor=(1.025, 1), loc='upper left')
plt.tight_layout()  # Fit everything inside the figure
plt.rcParams.update({'figure.autolayout': True})
plt.show()

print('Please proceed to the highground using the route displayed on the map. '
      'The route is ' + str(path_length) + ' km long.')
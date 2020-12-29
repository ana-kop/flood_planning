import pandas
import json
import networkx as nx
import rasterio.mask
from rasterio.plot import plotting_extent
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point, MultiPoint, box, LineString
import shapely.ops
import numpy as np

"""Import data"""
# Elevation data
path = "/Users/ak/Desktop/a_2/Material/"
background = rasterio.open(path + "/background/raster-50k_2724246.tif")
elevation = rasterio.open(path + "/elevation/SZ.asc")
# count = 1 // CRS = BNG
# rasterio.plot.show(elevation)

# ITN data
itn = pandas.read_json(path + "itn/solent_itn.json")
# print(itn['roadlinks'].head())
# print(list(itn))

# Nodes and links data
nodes = gpd.read_file(path + "roads/nodes.shp")
# print(nodes.head())
# print(list(nodes))
# with pandas.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified
#      print(nodes)

links = gpd.read_file(path + "roads/links.shp")

"""  Task 1:
Ask the user to input their current location as a British National Grid coordinate (easting and northing),
test whether the user is within a box (430000, 80000) and (465000, 95000); quit program if not
"""

# print("This program will help you find the the quickest route to walk to the highest point of land "
#      "within a 5km radius.")
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
# # test whether the user is within a box (430000, 80000) and (465000, 95000).
# inside = True
#
# # code for testing that user is within the (430000, 80000) and (465000, 95000)
#
# if inside is False:
#    print("Unable to assist in finding highest point of land.")

x = 433600  # 449923.625
y = 86600  # 89243.008

#bottom = 450000, 76000

# from Yolanda
# test whether the user is within a box (430000, 80000) and (465000, 95000).
# if x <= 430000 or x >= 465000 or y <= 80000 or y >= 95000:
#     print("Unable to assist in finding highest point of land because it's not within a box (430000, 80000) and "
#           "(465000, 95000).")
#     exit()
# else:
#     print("Thanks for inputting")

user_location = Point(x, y)  # Example user location
buf = user_location.buffer(5000)  # create 5km buffer polygon.

"""Task 2: Find highest point within 5km radius
"""
# https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html

buffer_gdf = gpd.GeoDataFrame({'geometry': buf}, index=[0], crs="EPSG:27700")  # put buffer polygon into geodataframe (gdf)

elevation_output_image, output_transform = rasterio.mask.mask(elevation,buffer_gdf['geometry'],\
                                                              nodata=-100,crop=True, filled=True)
out_meta = elevation.meta       # match the output image metadata with elevation metadata
out_meta.update({"driver": "GTiff",     # update metadata for elevation output img
                 "height": elevation_output_image.shape[1],
                 "width": elevation_output_image.shape[2],
                 "transform": output_transform})
# elevation.close()

with rasterio.open(path+"/elevation_output.tif", "w", **out_meta) as dest:   # write file
    dest.write(elevation_output_image)

radius = rasterio.open(path+"/elevation_output.tif", "r")       # 5km radius image file opened in read mode
radius_array= radius.read(1)                                    # radius read as numpy array

max_height = np.max(radius_array)   # find max value within buffer
pix_location_y, pix_location_x = np.where(radius_array == max_height)   # finds row/columns of pixel value
highest_points = []

for i in range (len(pix_location_x)):
    bng_pixel_location = radius.transform*(pix_location_x[i],pix_location_y[i]) # transforms (row,column) into (x,y)
    high_point = Point(bng_pixel_location) # create shapely point at highest point
    highest_points.append (high_point)    # append to list

"""Task 3.
Identify the nearest ITN node to the user and the nearest ITN node to the highest point identified in Task 2.
"""

# https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
# For every node in the nodes gdf, take the corresponding point from the geometry column and append it to nodes_list.
nodes_list = []
for i in range(len(nodes)):
    node_point = nodes.at[i, 'geometry']
    nodes_list.append(node_point)
nodes_locations = MultiPoint(nodes_list)  # pass the list of node points to a Multipoint shapely object


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

# Creating an empty multidigraph object
G = nx.MultiDiGraph()

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
        el = list(elevation.sample([(x, y)]))[0][0]
        road_coords_elev.append(el)

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
    time_to_cover = ((3 * length) / 250) + (elev / 10)

    # Adding edges with their corresponding weights
    G.add_edge(road_links[link_fid]['start'], road_links[link_fid]['end'], fid=link_fid, weight=time_to_cover)

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

    time_to_cover = ((3 * length) / 250) + (elev / 10)
    G.add_edge(road_links[link_fid]['end'], road_links[link_fid]['start'], fid=link_fid, weight=time_to_cover)

# Finding the shortest path using the dijkstra algorithm from networkx.
the_path = nx.dijkstra_path(G, nearest_node_user, nearest_node_highground, weight="weight")

# Converting the shortest path from FID codes to a linestring of coordinates
with open(path + "itn/solent_itn.json") as file:
    roadnodes = json.load(file)['roadnodes']

coords_path = []  # list of road points
for i in range(len(the_path)):
    osgb_code = the_path[i]
    coords = roadnodes[osgb_code]['coords']
    coords = tuple(coords)
    coords_path.append(coords)
shortest_path = LineString(coords_path)


"""Task 5: Output
"""
# bounding box dimensions (for plotting the background file):
minx = user_location.x - 5000  # 5km in each direction
maxx = user_location.x + 5000
miny = user_location.y - 5000
maxy = user_location.y + 5000
bbox = box(minx, miny, maxx, maxy)

bbox_gdf = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs="EPSG:27700") # puts bounding box into Geodataframe

# crop the background image to 10km x 10km using bounding box:
cropped_background_img, background_transform = rasterio.mask.mask(background,bbox_gdf['geometry'],crop=True)
background_meta = background.meta       # match cropped background img metadata with original background metadata

background_meta.update({"driver": "GTiff",           # update metadata for output img
                 "height": cropped_background_img.shape[1],
                 "width": cropped_background_img.shape[2],
                 "transform": background_transform})

with rasterio.open(path+"/background_output.tif", "w", **background_meta) as dest:  # write file
    dest.write(cropped_background_img)

background.close()

cropped_background = rasterio.open(path+"/background_output.tif", "r")   # opens cropped background
background_extent = plotting_extent(cropped_background) # extent used for plotting


# MASKING THE ELEVATION FILE TO PLOT ONLY THE CIRCULAR BUFFER
circle = np.ma.masked_where(radius_array == -100, radius_array)

# Masking background sea level values to re-colour
removed_sea = np.ma.masked_where(cropped_background.read(1)==65, cropped_background.read(1))

# PLOT EVERYTHING:
fig, ax = plt.subplots(figsize=(8,6))
ax.set_facecolor("paleturquoise")       # sets sea colour
ax.imshow(removed_sea, cmap="terrain", extent=background_extent) # plots cropped background
elevation_radius = ax.imshow(circle, alpha=0.55, aspect=1, extent=background_extent, cmap="viridis") #plots elevation
ax.plot(user_location.x, user_location.y, 'o', color='tomato', markersize=8, label="User Location")    # plot user location
ax.plot(highest_points[0].x,highest_points[0].y,'^',color='black',markersize=8,label="Highest Point") # plot highest point
cbar =plt.colorbar(elevation_radius, orientation='vertical', fraction=0.025, pad=0.12) #plot colour bar
cbar.set_label('Elevation (m)') #plot colour bar label
plt.arrow((minx+500), (maxy-1200),0,950, width=120, head_length=300, length_includes_head=True, facecolor="black", edgecolor= "black")   #plot north arrow
ax.plot([(maxx-2500), (maxx-500)], [(miny+475),(miny+475)], color="black", linewidth=1.5) # plot 2km scale bar:
plt.text(maxx-1750, miny+125, "2km", fontsize=8, fontweight='bold')   #plot scale bar label
plt.legend(bbox_to_anchor=(1.025, 1), loc='upper left') # plot legend

x, y = shortest_path.xy
ax.plot(x, y, color='black', zorder=1, label="shortest path")

plt.tight_layout() # fit everything inside the figure
plt.rcParams.update({'figure.autolayout': True})
plt.show() #show figure

print('Please proceed to highground at '
      + str(nearest_node_highground) + '. The closest road for you to go there starts at point '
      + str(nearest_node_user) + '.')
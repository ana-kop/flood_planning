import rasterio.plot
import pandas
import geopandas as gpd
from shapely.geometry import Point, MultiPoint
from shapely.ops import nearest_points

# Importing data

# Elevation data
path = "/Users/ak/Desktop/a_2/Material/"
background = rasterio.open(path + "/background/raster-50k_2724246.tif")
elevation = rasterio.open(path + "/elevation/SZ.asc")
# count = 1 // CRS = BNG
# rasterio.plot.show(elevation)

# ITN data
itn = pandas.read_json(path + "itn/solent_itn.json")
# print(itn.head())
# print(list(itn))

# Nodes and links data
nodes = gpd.read_file(path + "roads/nodes.shp")
# print(nodes.head())
links = gpd.read_file(path + "roads/links.shp")
# print(links.head())

'''  Task 1:
Ask the user to input their current location as a British National Grid coordinate (easting and northing), 
test whether the user is within a box (430000, 80000) and (465000, 95000); quit program if not
'''

print("This program will help you find the the quickest route to walk to the highest point of land "
     "within a 5km radius.")

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

''' Task 3:
Identify the nearest ITN node to the user and the nearest ITN node to the highest point identified in Task 2.
'''

x = 450000
y = 93500
location_user = Point(x, y)

location_highground = Point(448860, 92000)

# For every line in nodes, take point from geometry, append point to list, pass it to Multipoint
# https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
nodes_list = []
for i in range(len(nodes)):
    point = nodes.at[i, 'geometry']
    nodes_list.append(point)
node_locations = MultiPoint(nodes_list)

nearest_to_user = nearest_points(location_user, node_locations)
# print(nearest_geoms[0]) - original point
nearest_node_user = nearest_to_user[1]
#print(nearest_node_user)

nearest_to_highground = nearest_points(location_highground, node_locations)
nearest_node_highground = nearest_to_highground[1]
#print(nearest_node_highground)

''' Task 5:
Output
'''

print('Please proceed to highground at '
      + str(nearest_node_highground) + '. The closest road for you to go there starts at point '
      + str(nearest_node_user) + '.')

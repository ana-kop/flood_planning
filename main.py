import rasterio.plot
import pandas

# Importing data

# Elevation data
path = "/Users/ak/Desktop/a_2/Material/"
background = rasterio.open(path + "/background/raster-50k_2724246.tif")
elevation = rasterio.open(path + "/elevation/SZ.asc")
# count = 1 // CRS = BNG
# rasterio.plot.show(elevation)

# ITN data
itn = pandas.read_json(path + "itn/solent_itn.json")
print(itn.head())


# Task 1
# ask the user to input their current location as a British National Grid coordinate (easting and northing)

print("--This program will help you find the the quickest route to walk to the highest point of land "
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
inside = False

# code for testing that user is within the (430000, 80000) and (465000, 95000)

if inside is False:
    print("Unable to assist in finding highest point of land.")

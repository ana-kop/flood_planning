# ask the user to input their current location as a British National Grid coordinate (easting and northing)

print("--This program will help you find the the quickest route to walk to the highest point of land "
      "within a 5km radius.")
#  test
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

# code for testing that user is within the box

if inside is False:
    print("Unable to assist in finding highest point of land.")

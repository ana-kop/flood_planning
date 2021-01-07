'''  Task 1:
Ask the user to input their current location as a British National Grid coordinate (easting and northing), 
test whether the user is within a box (430000, 80000) and (465000, 95000); quit program if not
'''

print("This program will help you find the the quickest route to walk to the highest point of land "
     "within a 5km radius.")

# Use the while loop to guarantee that the user inputs Easting coordinate(x) as a number
while True:
    try:
        x = float(input("Please enter the Easting coordinate of your location: "))
        break
    # If user entered not a number and there is ValueError, then:
    except ValueError:
        print("Please input X-coordinate of your point again as a NUMBER.")

# Use the while loop to guarantee that the user inputs Northing coordinate(y) as a number
while True:
    try:
        y = float(input("Please enter the Northing coordinate of your location: "))
        break
    # If user entered not a number and there is ValueError, then:
    except ValueError:
        print("Please input Y-coordinate of your point again as a NUMBER.")

# test whether the user is within a box (430000, 80000) and (465000, 95000).
if x <= 430000 or x >= 465000 or y <= 80000 or y >= 95000:
    print("Unable to assist in finding highest point of land"
           "because it's not within a box (430000, 80000) and (465000, 95000).")
    exit()
else:
    print("Thanks for inputting")
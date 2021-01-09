import geopandas as gpd
import json
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import rasterio.mask
from rasterio.plot import plotting_extent
from shapely.geometry import Point, MultiPoint, box, LineString, MultiLineString
import shapely.ops
from shapely.ops import linemerge
from pyproj import Transformer
import what3words


def main(background_f, elevation_f, nodes_f, itn_f, shape_f, island_checker=False, allow_near_edge=False,
         input_options=False):
    """Definition of the main() function for the program"""

    """Task 1: User Input"""

    if input_options is False:  # The basic version of accepting user input
        print("This program will help you find the the quickest route to walk to the highest point of land "
              "within a 5km radius from your current location.")

        # Ask user to input their location and make a Point object with the coordinates
        try:
            x = float(input("Please enter the Easting coordinate of your location: "))
        # If user entered not a number and there is ValueError, then:
        except ValueError:
            x = float(input("Please input a NUMBER. Please enter the X-coordinate of your point: "))

        try:
            y = float(input("Please enter the Northing coordinate of your location: "))
        except ValueError:
            y = float(input("Please input a NUMBER. Please enter the Y-coordinate of your point: "))

        user_location = Point(x, y)

    # CREATIVITY TASK EXTENSION: giving user more options for inputting their location
    if input_options is True:
        # Ask user how they would like to input their location
        option = float(input('How would you like to enter your coordinates? '
                             ' \nFor British National Grid, enter 1; for latitude and longitude, enter 2;'
                             ' for what3words, enter 3: '))

        if option == 1:
            try:
                x = float(input("Please enter the Easting coordinate of your location: "))
            except ValueError:
                x = float(input("Please input a NUMBER. Please enter the X-coordinate of your point: "))

            try:
                y = float(input("Please enter the Northing coordinate of your location: "))
            except ValueError:
                y = float(input("Please input a NUMBER. Please enter the Y-coordinate of your point: "))

            user_location = Point(x, y)

        transformer = Transformer.from_crs("EPSG:4326", "EPSG:27700", always_xy=True)

        if option == 2:
            try:
                lat_x = float(input("Please enter the latitude of your location: "))
            except ValueError:
                lat_x = float(input("Please input a NUMBER. Please enter the X-coordinate of your point: "))

            try:
                lon_y = float(input("Please enter the longitude of your location: "))
            except ValueError:
                lon_y = float(input("Please input a NUMBER. Please enter the Y-coordinate of your point: "))

            # Transform the lat/lon coordinates into BNG using Transformer from pyproj
            output = transformer.transform(lat_x, lon_y)
            user_location = Point(output)

        if option == 3:
            word_1 = str(input("Please enter word 1: "))
            word_2 = str(input("Please enter word 2: "))
            word_3 = str(input("Please enter word 3: "))
            words = word_1 + '.' + word_2 + '.' + word_3

            geocoder = what3words.Geocoder("JWCSJO0N")  # Inputting the API key
            # Transform the 3 words into  lat/lon coordinates using the what3words API
            # https://developer.what3words.com/tutorial/python
            res = geocoder.convert_to_coordinates(words)
            x = res['square']['northeast']['lng']
            y = res['square']['northeast']['lat']
            output = transformer.transform(x, y)  # Transform the lat/lon coordinates into BNG using Transformer
            user_location = Point(output)

    # Test whether the user is within a box (430000, 80000) and (465000, 95000)
    if allow_near_edge is False:
        if user_location.x <= 430000 or user_location.x >= 465000 \
                or user_location.y <= 80000 or user_location.y >= 95000:
            print("Unable to assist in finding highest point of land.")
            exit()
        else:
            print("Input accepted, finding the shortest route to high ground...")

    if allow_near_edge is True:
        print("Input accepted, finding the shortest route to high ground..."
              " \n(Note: currently unable to display elevation map; only the shortest route will be displayed.)")

    # CREATIVITY TASK: testing whether user is outside the Isle of Wight
    if island_checker is True:
        island_gpd = gpd.read_file(shape_f)  # Read in the Isle of Wight shapefile
        island = island_gpd['geometry'][0]  # Get the MultiPolygon object with the boundary
        if island.contains(user_location) is False:
            print('Looks like you are not on the Isle of Wight. '
                  'This app is only intended to be used on the Isle of Wight.')
            exit()

    """Task 2: Highest Point Identification"""
    # Importing data
    background = rasterio.open(background_f)
    elevation = rasterio.open(elevation_f)

    # https://rasterio.readthedocs.io/en/latest/topics/masking-by-shapefile.html
    buf = user_location.buffer(5000)  # create 5km buffer polygon
    # Put buffer polygon into geodataframe
    buffer_gdf = gpd.GeoDataFrame({'geometry': buf}, index=[0], crs="EPSG:27700")

    # Mask elevation with buffer polygon and crop; points outside of the buffer set to value of -100
    elevation_output_image, output_transform = rasterio.mask.mask(elevation, buffer_gdf['geometry'],
                                                                  nodata=-100, crop=True, filled=True)

    out_meta = elevation.meta  # Match the output image metadata with elevation metadata
    out_meta.update({"driver": "GTiff",  # Update metadata for elevation output image
                     "height": elevation_output_image.shape[1],
                     "width": elevation_output_image.shape[2],
                     "transform": output_transform})

    with rasterio.open("Material/elevation_output.tif", "w", **out_meta) as dest:  # Write file
        dest.write(elevation_output_image)

    radius = rasterio.open("Material/elevation_output.tif", "r")  # 5km radius image file opened in read mode
    radius_array = radius.read(1)  # Read radius as numpy array

    max_height = np.max(radius_array)  # Find maximum value within buffer
    pix_location_y, pix_location_x = np.where(radius_array == max_height)  # Finds row/columns of pixel value
    highest_points = []

    def test_transformation(point):
        """Function to check that transformed coordinates fall within the bounds of the elevation file"""
        if point.x < 42500 or point.x > 470000 or point.y < 75000 or point.y > 100000:
            return True

    for i in range(len(pix_location_x)):
        # Transforms (row,column) into (x,y)
        bng_pixel_location = radius.transform * (pix_location_x[i], pix_location_y[i])
        high_point = Point(bng_pixel_location)  # Create shapely point at highest point
        if test_transformation(high_point):
            raise Exception(
                "Transformation Error: The determined coordinates do not fall within the bounds of the input "
                "file. Please ensure all files use the British National Grid CRS. ")
        else:
            highest_points.append(high_point)

    """Task 3: Nearest Integrated Transport Network"""

    # Read in the nodes shapefile into a gdf
    nodes = gpd.read_file(nodes_f)

    # For every node in the nodes gdf, check if it is less than 5 km away from the user_location;
    # append it to the nodes_list if it is.
    # https://automating-gis-processes.github.io/2017/lessons/L3/nearest-neighbour.html
    nodes_list = []
    for i in range(len(nodes)):
        node_i = nodes.at[i, 'geometry']
        if node_i.distance(user_location) < 5000:
            nodes_list.append(node_i)

    # Make a shapely MultiPoint object from the points in the nodes_list
    nodes_multipoint = MultiPoint(nodes_list)

    def nearest_itn(point, multipoint):
        """Definition of the function that takes location of a point of interest and a multipoint object containing
        the nodes, and returns the location of the ITN closest to the point of interest
        """
        nearest_nodes_tup = shapely.ops.nearest_points(point, multipoint)
        # Select the first ITN node from the tuple returned by nearest_points method
        # (nearest_nodes_tup[0] - original point)
        nearest_node = nearest_nodes_tup[1]
        return nearest_node

    nearest_node_user_coord = nearest_itn(user_location, nodes_multipoint)
    nearest_node_highground_coord = nearest_itn(highest_points[0], nodes_multipoint)

    # Finding the FID of a given ITN node: make a dictionary that allows mapping between the coordinate
    # location of an ITN and its FID
    itn_fid_dict = {}
    for i in nodes.index:
        geom = nodes['geometry'][i]
        fid_ = nodes['fid'][i]
        geom = str(geom.x) + ',' + str(geom.y)  # The point object is stored in the dictionary in a string format
        itn_fid_dict[geom] = fid_

    def coord_to_fid(itn_coord):
        """Definition of the function that takes coordinates of a given ITN node and returns its FID"""
        point_to_convert = str(itn_coord.x) + ',' + str(itn_coord.y)
        res_fid = itn_fid_dict[point_to_convert]
        return res_fid

    nearest_node_user = coord_to_fid(nearest_node_user_coord)
    nearest_node_highground = coord_to_fid(nearest_node_highground_coord)

    """Task 4: Shortest Path"""

    # Create an empty MultiDiGraph structure
    G = nx.MultiDiGraph()

    # Load the ITN data
    with open(itn_f) as f:
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
        for i in range(len(road_coords_elev) - 1):
            j = i + 1
            b = road_coords_elev[j]
            a = road_coords_elev[i]
            if b > a:
                elev = elev + (b - a)

        # Formula to calculate the weight of the node
        # - i.e. the time required to walk it according to the Naismith's rule.
        time_to_walk = ((3 * length) / 250) + (elev / 10)

        # Adding edges with their corresponding weights to the graph
        G.add_edge(road_links[link_fid]['start'], road_links[link_fid]['end'], fid=link_fid, weight=time_to_walk)

        # Calculate weights and add edges in the opposite direction
        road_coords_elev.reverse()
        for i in range(len(road_coords_elev) - 1):
            j = i + 1
            b = road_coords_elev[j]
            a = road_coords_elev[i]
            if b > a:
                elev = elev + (b - a)

        time_to_walk = ((3 * length) / 250) + (elev / 10)
        G.add_edge(road_links[link_fid]['end'], road_links[link_fid]['start'], fid=link_fid, weight=time_to_walk)

    # Finding the shortest path using the dijkstra algorithm from networkx.
    shortest_path = nx.dijkstra_path(G, nearest_node_user, nearest_node_highground, weight="weight")

    # Converting the FIDs returned by the algorithm into a linestring for plotting; getting path length.
    # https://gis.stackexchange.com/questions/223447/weld-individual-line-segments-into-one-linestring-using-shapely
    shortest_path_lines = []
    path_length = 0
    for i, c in enumerate(shortest_path[:-1]):
        fid = G.get_edge_data(c, shortest_path[i + 1])[0]['fid']
        coords = road_links[fid]['coords']
        # Convert this road segment to a LineString and append it to the list shortest_path_lines
        shortest_path_lines.append(LineString(coords))
        path_length = path_length + road_links[fid]['length']
    # Merge the LineStrings in the list into one LineString object
    shortest_path_final = linemerge(MultiLineString(shortest_path_lines))
    # Find the length of the shortest path
    path_length = round(path_length / 1000, 2)

    """Task 5: Map Plotting"""
    # Bounding box dimensions (for plotting the background file):
    minx = user_location.x - 5000  # 5km in each direction
    maxx = user_location.x + 5000
    miny = user_location.y - 5000
    maxy = user_location.y + 5000
    bbox = box(minx, miny, maxx, maxy)

    # Bounding box dimensions for plotting the 10kmx10km background file
    bbox_gdf = gpd.GeoDataFrame({'geometry': bbox}, index=[0], crs="EPSG:27700")

    # Crop the background image to 10km x 10km using bounding box:
    cropped_background_img, background_transform = rasterio.mask.mask(background, bbox_gdf['geometry'], crop=True)
    background_meta = background.meta  # Match cropped background img metadata with original background metadata

    background_meta.update({"driver": "GTiff",  # Update metadata for output img
                            "height": cropped_background_img.shape[1],
                            "width": cropped_background_img.shape[2],
                            "transform": background_transform})

    with rasterio.open("Material/background_output.tif", "w", **background_meta) as dest:  # write file
        dest.write(cropped_background_img)
    background.close()

    cropped_background = rasterio.open("Material/background_output.tif", "r")
    background_extent = plotting_extent(cropped_background)  # Background_extent used for plotting
    circle = np.ma.masked_where(radius_array == -100, radius_array)  # Masking elevation file to plot circular buffer

    # Masking background sea level values to re-colour
    removed_sea = np.ma.masked_where(cropped_background.read(1) == 65, cropped_background.read(1))

    # Plot everything:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor("paleturquoise")  # Set sea colour
    ax.imshow(removed_sea, cmap="terrain", extent=background_extent)  # Plot cropped background

    if allow_near_edge is False:  # Task 6 extension: plot elevation if inside the box
        elevation_radius = ax.imshow(circle, alpha=0.5, aspect=1, extent=background_extent, cmap="viridis")
        cbar = plt.colorbar(elevation_radius, orientation='vertical', fraction=0.025, pad=0.12)
        cbar.set_label('Elevation (m)')

    ax.plot(user_location.x, user_location.y, 'o', color='orangered', markersize=8, label="User Location")
    ax.plot(highest_points[0].x, highest_points[0].y, '^', color='black', markersize=8, label="Highest Point")
    plt.arrow((minx + 500), (maxy - 1200), 0, 950, width=120, head_length=300, length_includes_head=True,
              facecolor="black", edgecolor="black")  # Plot north arrow
    ax.plot([(maxx - 2500), (maxx - 500)], [(miny + 475), (miny + 475)], color="black",
            linewidth=1.5)  # Plot 2km scale bar
    plt.text(maxx - 1750, miny + 125, "2km", fontsize=8, fontweight='bold')  # Plot scale bar label

    # Plot the shortest path
    x, y = shortest_path_final.xy
    ax.plot(x, y, color='black', linewidth=2, zorder=1, label="shortest path")

    plt.legend(bbox_to_anchor=(1.025, 1), loc='upper left')
    plt.tight_layout()  # Fit everything inside the figure
    plt.rcParams.update({'figure.autolayout': True})
    plt.show()

    print('Shortest route found! Please proceed to the high ground using the route displayed on the map. '
          'The route is ' + str(path_length) + ' km long.')


if __name__ == "__main__":
    background_file = "Material/background/raster-50k_2724246.tif"
    elevation_file = "Material/elevation/SZ.asc"
    nodes_file = "Material/roads/nodes.shp"
    itn_file = "Material/itn/solent_itn.json"
    shape_file = "Material/shape/isle_of_wight.shp"

    main(background_file, elevation_file, nodes_file, itn_file, shape_file,
         island_checker=False, allow_near_edge=False, input_options=False)

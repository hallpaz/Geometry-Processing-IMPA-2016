import numpy as np
import triangle
import os


data_folder = "data"
images_folder = os.path.join("images", "crust")

def crust(points: list, draw_voronoi = False, draw_delaunay = False)-> list:
    points = np.array(points)
    ordered_points = []

    if draw_voronoi:
        draw_voronoi_diagram()

    if draw_delaunay:
        draw_delaunay_diagram()

    return ordered_points


def main():
    pass

if __name__ == '__main__':
    main()

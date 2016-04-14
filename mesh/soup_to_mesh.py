import argparse
import os

from DataStructures import Vertex, Triangle


data_folder = "data/"

def read_OFF(off_file):
    vertexBuffer = []
    indexBuffer = []
    with open(off_file, "r") as modelfile:
        first = modelfile.readline().strip()
        if first != "OFF":
            raise(Exception("not a valid OFF file ({})".format(first)))

        parameters = modelfile.readline().strip().split()
        min_distance = BIG
        minx = BIG
        miny = BIG
        minz = BIG
        maxx = -BIG
        maxy = -BIG
        maxz = -BIG

        if len(parameters) < 2:
            raise(Exception("OFF file has invalid number of parameters"))

        for i in range(int(parameters[0])):
            coordinates = modelfile.readline().split()
            X, Y, Z = float(coordinates[0]), float(coordinates[1]), float(coordinates[2])
            vertexBuffer.append(Vertex(X, Y, Z))
            if X < minx:
                minx = X
            if X > maxx:
                maxx = X
            if Y < miny:
                miny = Y
            if Y > maxy:
                maxy = Y
            if Z < minz:
                minz = Z
            if Z > maxz:
                maxz = Z

        for i in range(int(parameters[1])):
            indices = modelfile.readline().split()
            indexBuffer.append(Triangle(int(indices[1]), int(indices[2]), int(indices[3])))

        bounding_box = BoundingBox(minx, miny, minz-1, maxx, maxy, maxz+1)

    return vertexBuffer, indexBuffer, bounding_box

def soup_to_mesh(filename: str, output_file: str, dimension = 2):
    vertices = []
    indices = []
    with open(filename, 'r') as soup_file:
        soup_lines = soup_file.readlines()
        for line in soup_lines:
            coordinates = line.split()
            t = []
            for i in range(3):
                v = Vertex(float(coordinates[0+dimension*i]), float(coordinates[1+dimension*i]), 0 if dimension == 2 else float(coordinates[2+dimension*i]))
                if v in vertices:
                    # repeated vertex, reetrieve index
                    index = vertices.index(v)
                    t.append(index)
                else:
                    # new vertex, new index
                    t.append(len(vertices))
                    vertices.append(v)
            indices.append(Triangle(t[0], t[1], t[2]))


    vertices = [str(v) for v in vertices]
    indices = [str(i) for i in indices]

    with open(output_file, 'w') as meshfile:
        meshfile.write(
        '''OFF
        %d %d 0
        %s%s
        '''%(len(vertices),len(indices), "".join(vertices), "".join(indices)))


### Complete this ###
def share_edge(edge: set, triangle: list):

    if edge.intersection(triangle) == edge

def face_association_from_OFF(filename: str):
    vertices, indices = read_OFF(filename)


    indices = [set([t.a, t.b, t.c]) for t in indices]
    for t in indices:
        [j for j in indices if ]
        neighbors.append()
#####################

def main(input_file, dimension = 2):
    dot_index = input_file.rfind(".")
    output_file = (input_file[:dot_index] if dot_index > 0 else input_file.strip()) + ".off"
    soup_to_mesh(input_file, output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 2D points (x, y) ... The output is written to a txt file
    ''')
    parser.add_argument('input_file', help='input list of points with X Y per line (format: txt)')
    args = parser.parse_args()
    if args.input_file == "all":
        rec_all_data(data_folder)
    else:
        main(os.path.join(data_folder, args.input_file))

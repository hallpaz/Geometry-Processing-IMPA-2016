import numpy as np


# returns a list of vertices and a list of triangles (both represented as numpy arrays)
def read_OFF(off_file):
    vertexBuffer = []
    indexBuffer = []
    with open(off_file, "r") as modelfile:
        first = modelfile.readline().strip()
        if first != "OFF":
            raise(Exception("not a valid OFF file ({})".format(first)))

        parameters = modelfile.readline().strip().split()

        if len(parameters) < 2:
            raise(Exception("OFF file has invalid number of parameters"))

        for i in range(int(parameters[0])):
            coordinates = modelfile.readline().split()
            vertexBuffer.append([float(coordinates[0]), float(coordinates[1]), float(coordinates[2])])

        for i in range(int(parameters[1])):
            indices = modelfile.readline().split()
            indexBuffer.append([int(indices[1]), int(indices[2]), int(indices[3])])

    return np.array(vertexBuffer), np.array(indexBuffer)

# receives a list of vertices and a list of indices (both as numpy arrays)
def write_OFF(output_file, vertices, indices):
    '''vertices and indices are lists of strings'''

    # converts indices and vertices to a string representation
    str_vertices = ["{} {} {}\n".format(v[0], v[1], v[2]) for v in vertices]
    str_indices = ["3 {} {} {}\n".format(i[0], i[1], i[2]) for i in indices]
    with open(output_file, 'w') as meshfile:
        meshfile.write(
        '''OFF
        %d %d 0
        %s%s
        '''%(len(str_vertices),len(str_indices), "".join(str_vertices), "".join(str_indices)))


# returns a list of numpy arrays
def read_points(filename: str):
    points = []
    with open(filename) as myfile:
        file_lines = myfile.readlines()
        for line in file_lines:
            content = line.split()
            content = [float(n) for n in content]
            # each element is a numpy array
            points.append(content)
    return np.array(points)


#points is a list of numpy arrays
def write_points(points:list, filename:str):
    if len(points) == 0:
        return None
    with open(filename, "w") as myfile:
        for point in points:
            if len(point) == 2:
                myfile.write("{} {}\n".format(point[0], point[1]))
            elif len(point) == 3:
                myfile.write("{} {} {}\n".format(point[0], point[1], point[3]))
            else:
                raise Exception("Points should have dimension 2 or 3")

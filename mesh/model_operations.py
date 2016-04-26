from DataStructures import Vertex, Triangle

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
            X, Y, Z = float(coordinates[0]), float(coordinates[1]), float(coordinates[2])
            vertexBuffer.append(Vertex(X, Y, Z))

        for i in range(int(parameters[1])):
            indices = modelfile.readline().split()
            indexBuffer.append(Triangle(int(indices[1]), int(indices[2]), int(indices[3])))

    return vertexBuffer, indexBuffer


def write_OFF(output_file, vertices, indices):
    '''vertices and indices are lists of strings'''
    vertices = [str(v) for v in vertices]
    indices = [str(i) for i in indices]
    with open(output_file, 'w') as meshfile:
        meshfile.write(
        '''OFF
        %d %d 0
        %s%s
        '''%(len(vertices),len(indices), "".join(vertices), "".join(indices)))

from geometrypack.dataio import read_OFF, write_OFF
import numpy as np



vertices, indices = read_OFF("sphere2.off")
vertices = [v/(np.linalg.norm(v)) for v in vertices]
write_OFF("unitsphere2.off", vertices, indices)

import heapq
import math
from enum import Enum

class Color(Enum):
    White = 1
    Gray = 2
    Black = 3

class Node:
    def __init__(self, identifier = None):
        self.neighbors = []
        self.costs = []
        self.id = identifier
        self.parent = None
        self.color = Color.White
        self.weight = float("inf")

    def __str__(self):
        representation = "id " + str(self.id) + "\n"
        for i, j in zip(self.neighbors, self.costs):
            representation += "(" + str(i) + ", " + str(j) + ") "
        return representation

    def add(self, neighbor: int, cost = 1):
        self.neighbors.append(neighbor)
        self.costs.append(cost)

    def __lt__(self, other):
         return (self.weight < other.weight)
    def __le__(self, other):
         return (self.weight <= other.weight)
    def __eq__(self, other):
         return (self.weight == other.weight)
    def __ne__(self, other):
         return (self.weight != other.weight)
    def __gt__(self, other):
         return (self.weight > other.weight)
    def __ge__(self, other):
         return (self.weight >= other.weight)

#A min priority-queue (not necessary)
class Min_Priority_Queue():
    def __init__(self, array = []):
        self.storage = array
        self.build_heap(self.storage)

    def min_heapify(self, index:int):
        """Makes the ith element of a list to respect the min heap property"""

        min_index = index
        left = 2*index+1
        if((left < self.size()) and (self.storage[left] <= self.storage[min_index])):
            min_index = left
        right = 2*index+2
        if((right < self.size()) and (self.storage[right] <= self.storage[min_index])):
            min_index = right
        if(min_index != index):
            aux_bucket = self.storage[index]
            self.storage[index] = self.storage[min_index]
            self.storage[min_index] = aux_bucket

            self.min_heapify(min_index)
        return

    def decrease_key(self, index, key):
        pass
    def add(element):
        pass

    def validate_heap(self)->str:

        for i in range(heap_size//2):
            left = 2*i+1
            right = 2*i+2
            if(left < self.size() and self.storage[i] > self.storage[2*i+1]):
                print("HEAP PROPERTY VIOLATION", self.size(), self.storage)
                break
            if(right < self.size() and self.storage[i] > self.storage[2*i+2]):
                print("HEAP PROPERTY VIOLATION", self.size(), self.storage)
                break
        return "HEAP IS OK"

    def build_heap(self, array:list):
        """Builds a heap of minimum from a list"""

        for i in range(len(array)//2, -1, -1):
            self.min_heapify(i)

    def size(self):
        return len(self.storage)
    def empty(self):
        return bool(self.storage)


class DepthMap:
    def __init__(self, points, width, height):
        self.points = points
        self.width = width
        self.height = height


class Vertex:
    def __init__(self, x, y, z, value = None):
        self.x = x
        self.y = y
        self.z = z
        if value is None:
            self.value = 1
        else:
            self.value = value

    def __add__(self, rhs):
        return Vertex(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)

    def invert(self):
        return Vertex(-self.x, -self.y, -self.z)

    def __sub__(self, rhs):
        return Vertex(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)

    def __mul__(self, rhs):
        return self.x * rhs.x + self.y * rhs.y + self.z * rhs.z

    def __truediv__(self, scalar):
        return Vertex(self.x/scalar, self.y/scalar, self.z/scalar)

    def scalar_mult(self, scalar):
        return Vertex(self.x*scalar, self.y*scalar, self.z*scalar)

    def normalize(self):
        norm = math.sqrt(self*self)
        if(norm < 0.000001):
            print("Trying to normalize 0 vector")
            return Vertex(0, 0, 0)
        return Vertex(self.x/norm, self.y/norm, self.z/norm)

    def distance(lhs, rhs):
        return math.sqrt((lhs.x - rhs.x)**2 + (lhs.y - rhs.y)**2 + (lhs.z - rhs.z)**2)

    def __eq__(self, rhs):
        if Vertex.distance(self, rhs) < 0.001:
            return True
        return False

    def cross(lhs, rhs):
        x = (lhs.y*rhs.z - lhs.z*rhs.y)
        y = (lhs.z*rhs.x - lhs.x*rhs.z)
        z = (lhs.x*rhs.y - lhs.y*rhs.x)
        return Vertex(x, y, z)

    def __str__(self):
        return "%f %f %f\n"%(self.x,self.y,self.z)

    def repr(self):
        return str(self)

class Triangle:
    #triangle indices
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

    def __str__(self):
        return "3 %d %d %d\n"%(self.a, self.b, self.c)

    def __repr__(self):
        return "{} {} {}".format(self.a, self.b, self.c)


class GridCell:
    def __init__(self, vertices):
        if(len(vertices) != 8):
            raise("Wrong number of vertices ({}) for cell".format(len(vertices)))
        self.vertices = vertices

    def __repr__(self):
        return "cell"

# ------------------------------------------------------------------
# Priority Queue implementation from Python docs
class MinPriorityQueue():
    REMOVED = '<removed-point>'      # placeholder for a removed point
    counter = 0     # unique sequence count

    def __init__(self):
        self.data = []
        self.entry_finder = {} # mapping of points to entries

    def add_point(self, point, priority=0):
        'Add a new point or update the priority of an existing point'
        if point in self.entry_finder:
            self.remove_point(point)
        count = MinPriorityQueue.counter
        MinPriorityQueue.counter += 1
        entry = [priority, count, point]
        self.entry_finder[point] = entry
        #print(self.data, entry)
        heapq.heappush(self.data, entry)

    def remove_point(self, point):
        'Mark an existing point as REMOVED.  Raise KeyError if not found.'
        entry = self.entry_finder.pop(point)
        entry[-1] = MinPriorityQueue.REMOVED

    def pop_point(self):
        'Remove and return the lowest priority point. Raise KeyError if empty.'
        while self.data:
            priority, count, point = heapq.heappop(self.data)
            if point is not MinPriorityQueue.REMOVED:
                del self.entry_finder[point]
                return (priority, point)
        raise KeyError('pop from an empty priority queue')

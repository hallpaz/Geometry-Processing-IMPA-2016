from Point2D import Point2D, CCW_test
import math


class Segment():

    def readFromFile(filename:str)->list:
        segments = []
        with open(filename) as myfile:
            file_lines = myfile.readlines()
            for line in file_lines[1:]:
                coordinates = line.split()
                if(len(coordinates) == 4):
                    p1 = Point2D(float(coordinates[0]), float(coordinates[1]))
                    p2 = Point2D(float(coordinates[2]), float(coordinates[3]))
                    segments.append(Segment(p1, p2))
                else:
                    print("Weird coordinates result is undefined\n")
        return segments

    def __init__(self, pointA:Point2D, pointB:Point2D, name = None):
        #first is always the leftmost
        if pointA.x < pointB.x:
            self.first = pointA
            self.second = pointB
        else:
            self.first = pointB
            self.second = pointA

        if name is not None:
            self.name = name
            return
        if pointA.name and pointB.name:
            self.name = pointA.name + pointB.name
        else:
            self.name = None

    def __str__(self):
        representation = "{0}, {1}".format(str(self.first), str(self.second))
        if self.name is not None:
            representation += " " + self.name
        return representation

    def intersects(self, other)->bool:
        d1 = CCW_test(other.second, self.first, other.first)
        d2 = CCW_test(other.second, self.second, other.first)
        d3 = CCW_test(self.second, other.first, self.first)
        d4 = CCW_test(self.second, other.second, self.first)

        if( ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0))):
            return True
        elif d1 == 0 and other.contains(self.first):
            return True
        elif d2 == 0 and other.contains(self.second):
            return True
        elif d3 == 0 and self.contains(other.first):
            return True
        elif d4 == 0 and self.contains(other.second):
            return True

        return False

    def ray_intersection(self, other)->float:
        x1, y1 = self.first.x, self.first.y
        x2, y2 = self.second.x, self.second.y
        x3, y3 = other.first.x, other.first.y
        x4, y4 = other.second.x, other.second.y

        d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
        if abs(d) < 0.00001:
            return None

        xi = ((x3-x4)*(x1*y2-y1*x2)-(x1-x2)*(x3*y4-y3*x4))/d
        yi = ((y3-y4)*(x1*y2-y1*x2)-(y1-y2)*(x3*y4-y3*x4))/d

        return Point2D(xi,yi)

        return None

    def contains(self, point)->bool:
        if (min(self.first.x, self.second.x) <= point.x and
        point.x <= max(self.first.x, self.second.x) and
        min(self.first.y, self.second.y) <= point.y and
        point.y <= max(self.first.x, self.second.y)):
            return True
        return False

    def distance_to_point(self, point)->float:
        length_squared = Point2D.dot_product(self.first - self.second, self.first - self.second)
        if length_squared < 0.00000001:
            return Point2D.distance(self.first, point)

        u = ((point.x - self.first.x)*(self.second.x - self.first.x) + (point.y - self.first.y)*(self.second.y - self.first.y))/length_squared
        projection = self.first + u*(self.second - self.first)
        return Point2D.distance(point, projection)


    def __lt__(self, other):
        """This method checks whether the current segment is located below the other segment or not.
        We assume that the segments don't intersect each other, because if they do the intersection
        happens after the given x (otherwise the algorithm would have returned) """

        if CCW_test(self.second, other.first, self.first) > 0:
            return True
        else:
            return False

    def isBelow(self, other, x:float)->bool:
        """This method checks whether the current segment is located below the other segment or not.
        We assume that the segments don't intersect each other, because if they do the intersection
        happens after the given x (otherwise the algorithm would have returned) """

        if CCW_test(self.second, other.first, self.first) > 0:
            return True
        else:
            return False


if __name__ == "__main__":
    print("Module called as main")

    p1 = Point2D(0, 0)
    p2 = Point2D(1, 1)
    p3 = Point2D(1, 0)
    p4 = Point2D(0, 1)

    segA = Segment(p1, p2)
    segB = Segment(Point2D(0.5, 0.5), p4)

    print( segA.intersects(segB) )

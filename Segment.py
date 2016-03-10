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

        t = max(0, min(1, Point2D.dot_product(point - self.first, self.second - self.first)))
        projection = self.first + t*(self.second - self.first)
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

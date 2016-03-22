from Point2D import Point2D
from enum import IntEnum


class Direction(IntEnum):
    """docstring for Direction"""
    East = 0
    South = 1
    West = 2
    North = 3


def pavlids_boundary_tracing(data):
    start_col = 0
    start_row = 0
    found = False
    seen = 0
    for row in range(len(data)):
        for col in range(len(data[0])):
            seen += 1
            if data[row][col]:
                start_col = col
                start_row = row
                found = True
                print(row, col, data[row][col], len(data), len(data[0]))
                break
        if found:
            break

    if not found:
        return None
    current_direction = Direction.East
    boundary = [Point2D(start_row, start_col)]

    def p1(direction, i, j):
        if direction is Direction.North:
            return i-1, j-1
        elif direction is Direction.West:
            return i+1, j-1
        elif direction is Direction.South:
            return i+1, j+1
        elif direction is Direction.East:
            return i-1, j+1
        else:
            raise "Invalid Direction"

    def p2(direction, i, j):
        if direction is Direction.North:
            return i-1, j
        elif direction is Direction.West:
            return i, j-1
        elif direction is Direction.South:
            return i+1, j
        elif direction is Direction.East:
            return i, j+1
        else:
            raise "Invalid Direction"

    def p3(direction, i, j):
        if direction is Direction.North:
            return i-1, j+1
        elif direction is Direction.West:
            return i-1, j-1
        elif direction is Direction.South:
            return i+1, j-1
        elif direction is Direction.East:
            return i+1, j+1
        else:
            raise "Invalid Direction"

    rotations = 0
    row, col = start_row, start_col
    while rotations <= 3:
        #print("it", row, col)
        i, j = p1(current_direction, row, col)
        try:
            next_point = data[i][j]
        except Exception as e:
            print("row: ", row, "col: ", col)
            print(i, j, "direction: ", current_direction, "error: ", e)
            raise
        if next_point:
            if i == start_row and j == start_col:
                break
            boundary.append(Point2D(i*10, j*10))
            row, col = i, j
            current_direction = Direction((current_direction - 1)%4)
            rotations = 0
        else:
            i, j = p2(current_direction, row, col)
            try:
                next_point = data[i][j]
            except Exception as e:
                print("row: ", row, "col: ", col)
                print(i, j, "direction: ", current_direction, "error: ", e)
                raise
            if next_point:
                if i == start_row and j == start_col:
                    break
                boundary.append(Point2D(i*10, j*10))
                row, col = i, j
                rotations = 0
                #current_direction = Direction((current_direction)%4)
            else:
                i, j = p3(current_direction, row, col)

                try:
                    next_point = data[i][j]
                except Exception as e:
                    print("row: ", row, "col: ", col)
                    print(i, j, "direction: ", current_direction, "error: ", e)
                    raise
                if next_point:
                    if i == start_row and j == start_col:
                        break
                    boundary.append(Point2D(i*10, j*10))
                    row, col = i, j
                    rotations = 0
                else:
                    current_direction = Direction((current_direction+1)%4)
                    rotations += 1
    print("rotations: ", rotations)
    print("row: ", row, "col: ", col)
    return boundary

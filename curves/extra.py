def get_samples_from_image(image, num_samples = 100):
    minx = BIG
    maxx = -BIG
    miny = BIG
    maxy = -BIG
    samples = []
    for j in range(image.height):
        for i in range(image.width):
            if image.getpixel((i,j)) == 0: #is black
                if i < minx:
                    minx = i
                if i > maxx:
                    maxx = i
                if j < miny:
                    miny = j
                if j > maxy:
                    maxy = j

    points = [Point2D(minx, miny), Point2D(minx, maxy), Point2D(maxx, maxy), Point2D(maxx, miny)]
    draw_segments(points, "original.eps")

    maxx += 1
    minx -= 1
    maxy += 1
    miny -= 1
    print((minx, miny), (maxx, maxy))
    center = Point2D(minx + (maxx - minx)/2, miny + (maxy - miny)/2)
    radius = math.sqrt(((maxx - minx)/2)**2 + ((maxy - miny)/2)**2)
    print(center, radius)
    inc = 2*math.pi/num_samples

    seg_list = []

    for i in range(num_samples):
        pos = center + radius*Point2D(math.cos(i*inc), math.sin(i*inc))
        sample = getsample(image, pos, center)
        seg_list.append(Segment(pos, center))
        if sample is not None:
            samples.append(sample)
            print("sample: ", i)

    draw_segment_set(seg_list, "segmentos.eps")
    return samples

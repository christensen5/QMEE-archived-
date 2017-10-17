import csv
import sys
import math
import scipy

def computeHeight(degrees, distance):
    radians = math.radians(degrees)
    height = distance * math.tan(radians)
    print("Tree height is:" + str(height))

    return height

def main(argv):
    # read i
    input_filename = '../Data/' + str(argv[1])
    f = open(input_filename, 'rb')
    namestr = str(argv[1])
    g = open('../Results/' + namestr[:-4] + '_treeheights_python.csv', 'wb')
    trees = csv.reader(f)
    csvwrite = csv.writer(g)
    next(trees) # skip header line

    csvwrite.writerow(['', 'Species', 'Distance.m', 'Angle.degrees', 'Tree.Height.m'])
    i = 1
    for row in trees:
        csvwrite.writerow([i, row[0], row[1], row[2], computeHeight(float(row[2]), float(row[1]))])
        i += 1

    return 0

if __name__ == "__main__":
    main(sys.argv)
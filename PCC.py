# -*- coding: utf-8 -*-

#created by Hunter Whitlock, last edited on 4/11/2025

from Classes import *
from Imports import *
from Settings import *

def importData(data_path) -> list:
    output = []

    pathList = Path(data_path).glob('**/*.txt')
    for path in pathList: #for each file that ends in .txt in the data_path directory (relative file path)
        print("reading: " + str(path))
        with open(path, 'r') as file: #opens the file for reading

            output.append(copy.deepcopy(FileData(filename=str(path).split("/")[-1]))) #creates new filedata item with name of the file currently being read

            data = pd.read_csv(str(path), header=None).values.tolist()

            for i in range(len(data)): #splits each string of XY coords into a list containing an X and Y coord
                data[i] = data[i][0].rsplit("\t")
                data[i] = (float(data[i][0]), float(data[i][1]))

            #necessary to prevent file concatenation
            output[-1].XY.clear()

            for i in range(len(data)): #creates a deep copy that is appended to the most recently created FileData
                output[-1].XY = copy.deepcopy(data)
            file.close()

    return output

def writeDataToCSV(dataList, outputDirectory):

    csvData = []

    #for each file read, adds values to csvData to be written to .csv
    for i in range(len(dataList)):
        rowData = {"Filename" : dataList[i].filename}

        rowData["Curvature"] = dataList[i].curvature
        rowData["Best Fit Type"] = dataList[i].bestFitType
        rowData["Parabolic MAPE"] = dataList[i].parabola.meanAbsolutePercentageError
        rowData["Elliptic MAPE"] = dataList[i].ellipse.meanAbsolutePercentageError
        csvData.append(copy.deepcopy(rowData))
    
    #writes data to csv
    with open(outputDirectory + "TotalCurvatures.csv", "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[item for item in csvData[0].keys()])
        writer.writeheader()
        writer.writerows(csvData)
        csvfile.close()


def main():
    dataList = importData(data_path=DATA_DIRECTORY)

    playingWith = Bezier()

    # controlPoints = [(83,140), (94,117), (117,97), (148,90), (173,97), (191,118), (199,140)]
    # numControlPoints = 7
    # weights = [1, 1, 1, 1, 1, 1, 1]


    #out1, out2 = playingWith.rationalBezierExpression(numControlPoints, controlPoints, weights, do_division=True)
    #print(str(sp.latex(out1)) + "\n" + str(str(sp.latex(out2))))

    playingWith.fitCurve(dataList[0])

    for i in range(len(dataList)):
        dataList[i].findCurvature() #iterates through list, calculating curvature
        print(dataList[i].filename)
        # dataList[i].printXY()
        dataList[i].parabola.printFormula()
        #dataList[i].ellipse.printFormula()

    dataList[0].parabola.printRotatedXY()


    writeDataToCSV(dataList,DATA_DIRECTORY)
    print("Curvatures stored in the .csv in the data Directory")
    
    return 0
main()
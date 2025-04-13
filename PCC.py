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

            data = pd.read_csv(str(path)).values.tolist()

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

#WIP
def curvatureTests():

    def circle(x,r):
        return np.sqrt(r**2-x**2)

    data = []
    model = None


    #generates (x,y) tuples along the upper half of the unit circle
    for num in np.linspace(-1,1, num=10):
        data.append((num,circle(num,1)))








    return 1

def main():
    dataList = importData(data_path=DATA_DIRECTORY)

    for i in range(len(dataList)):
        dataList[i].findCurvature() #iterates through list, calculating curvature
        print(dataList[i].filename)
        # dataList[i].printXY()
        # dataList[i].parabola.printFormula()
        dataList[i].ellipse.printFormula()

    # print(str(sp.sympify(dataList[3].ellipse.curvatureFormula[0]).subs(dataList[3].ellipse.coefficients).subs("t",2).evalf()) + "," + str(
    #     sp.sympify(dataList[3].ellipse.curvatureFormula[1]).subs(dataList[3].ellipse.coefficients).subs("t", 2).evalf()))

    writeDataToCSV(dataList,DATA_DIRECTORY)
    print("Curvatures stored in the .csv in the data Directory")
    
    return 0

main()
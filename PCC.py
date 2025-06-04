from Classes import *
from Imports import *

DATA_DIRECTORY = "./Data/" #where the .txt's are located and where the .csv will be written

def importData(data_path) -> list:
    output = []

    pathList = Path(data_path).glob('**/*.txt')
    for path in pathList: #for each file that ends in .txt in the data_path directory (relative file path)
        print("Reading: " + str(path))
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
        rowData["x_curve"] = sp.latex(sp.S(dataList[i].bezier.curve["x_curve"]))
        rowData["y_curve"] = sp.latex(sp.S(dataList[i].bezier.curve["y_curve"]))
        csvData.append(copy.deepcopy(rowData))
    
    #writes data to csv
    with open(outputDirectory + "TotalCurvatures.csv", "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[item for item in csvData[0].keys()])
        writer.writeheader()
        writer.writerows(csvData)
        csvfile.close()


def main():
    # # testing degree elevation
    # t, u = sp.symbols("t, u", real=True)
    # test = Bezier()
    # interpolation_points = [(1,0), (12/13,5/13), (0,1)]
    # t_vals = [0,1/3,1]
    # weights = [2, 13/6, 1/2]
    #
    # x_curve, y_curve = test.barycentric_expression(interpolation_points, weights, t_vals)
    # print(sp.latex(sp.S(x_curve.replace('j','t'))) + "\n" + sp.latex(sp.S(y_curve.replace('j','t'))))
    #
    # x_curve, y_curve = sp.parse_expr(str(x_curve).replace('j', 't'), local_dict={'t': t}), sp.parse_expr(str(y_curve).replace('j', 't'),local_dict={'t': t})
    # dx_dt, dy_dt = x_curve.diff(t), y_curve.diff(t)
    # #arc_length = sp.integrate((dx_dt**2 + dy_dt**2)**0.5, t)
    # print(sp.latex(arc_length))

    test = Tests()
    test.All_Tests()


    dataList = importData(data_path=DATA_DIRECTORY)
    for i in range(len(dataList)):
        print("Analyzing: " + dataList[i].filename)
        dataList[i].findCurvature() #iterates through list, calculating curvature


    writeDataToCSV(dataList,DATA_DIRECTORY)
    print("Curvatures stored in the .csv in the data directory")
    
    return 0
main()
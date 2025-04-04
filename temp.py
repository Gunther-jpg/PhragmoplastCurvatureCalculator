# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#created by Hunter Whitlock on 3/17/2025 with code borrowed from Mark Dickinson

import math
import warnings
import copy
import csv
import scipy.optimize
import sympy as sp
import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass

warnings.filterwarnings('ignore')

DATA_DIRECTORY = "./Data/" #where the files are to be stored
global DATA_LIST; DATA_LIST = []
ORIGIN = (0,0)

@dataclass #stores data associated with each xy-coordinate file
class fileData:
    filename: str
    theta: float
    coefficients: list
    meanSquaredError: float
    X: list
    Y: list

    def __init__(self, filename="None", XDataList = [], YDataList = [], theta = 0, coefficients = [], meanSquaredError = -1):
        self.filename = filename
        self.X = XDataList
        self.Y = YDataList
        self.theta = theta
        self.coefficients = coefficients
        self.meanSquaredError = meanSquaredError
    
    def print(self):
        print("filename: " + self.filename)
        print("\ntheta: " + str(self.theta))
        print("\ncoefficients: " + str(self.coefficients[0]) + " " + str(self.coefficients[1]) + " " + str(self.coefficients[2]))
        print("\nmean squared error: " + str(self.meanSquaredError))
    
    def printXY(self):
        for i in range(len(self.X)):
            print(str(self.X[i]) + "\t" + str(self.Y[i]))

      
def importData(data_path) -> list:
    output = []

    pathList = Path(data_path).glob('**/*.txt')
    for path in pathList:
        print("reading: " + str(path))
        with open(path, 'r') as file:
            
            #currentData = fileData(str(path)) #creates new fileData object
            output.append(copy.deepcopy(fileData(str(path).split("/")[-1])))

            #sets data equal to a list of strings in the form of "$XPoint $YPoint"
            #data = [row[0] for row in csv.reader(file,delimiter=' ')] #https://stackoverflow.com/questions/17056526/what-is-the-easiest-way-to-read-the-text-file-delimited-by-tab-in-python
        
            data = pd.read_csv(str(path)).values.tolist()

            for i in range(len(data)): #splits each string of XY coords into a list containing an X and Y coord
                data[i] = data[i][0].rsplit("\t")
            
            #necessary to prevent file concatenation
            output[-1].X.clear()
            output[-1].Y.clear()

            for i in range(len(data)): #takes each point represented as a list in the form of ["$XCoord", "$YCoord"] and appends it as a float to the currentData class
                output[-1].X.append(float(data[i][0]))
                output[-1].Y.append(float(data[i][1]))
            #output[-1].printXY()
            #print("\n")

            file.close()

            #output.append(currentData) #appends the currentData to output and repeats the loop
            

    return output


########################################################################################
#    Author: Mark Dickinson - https://stackoverflow.com/users/270986/mark-dickinson
#    Date: 12/19/2015
#    Availability: https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
########################################################################################

def rotate(origin, point, angle): #https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy


def generalParabola(XYPair, a, b, c, d, e, f):
    x,y = XYPair[0], XYPair[1]

    output = a*(x**2) + b*x*y + c*y**2 + d*x + e*y + f #general equation for a parabola of any rotation/positon

    return output

def standardParabola(x, a, b, c):
    return a*x**2 + b*x + c

def ellipse(x, a, b, c, d, f):
    X,a0,a1,a2,a3,a4 = sp.symbols('X,a0,a1,a2,a3,a4') 
    equation = a0 * sp.cos(sp.asin((X-a1)/-a2) + a3) + a4
    equation = equation.subs({X: x, a0: a, a1: b, a2: c, a3: d, a4: f})
    
    return float(equation.evalf())
    
    

def rotateCurve(theta, XYData):
    theta = math.radians(theta)

    INDEX = DATA_LIST.index(XYData)
    meanSquaredError = 0
    predicted_values, rotatedX, rotatedY = [], [], []

    for i in range(len(XYData.X)): #rotates the xy data by a given theta, in radians
        x, y = XYData.X[i], XYData.Y[i]
        rotatedX.append(rotate(ORIGIN, (x,y), theta)[0])
        rotatedY.append(rotate(ORIGIN, (x,y), theta)[1])
        
        #print(str(rotatedX[i]) + "\t" + str(rotatedY[i]))
    #print("\n")
        
    results = scipy.optimize.curve_fit(f=ellipse, xdata=rotatedX, ydata=rotatedY)
    results = results[0]
    
    #calulates the mean squared error
    for i in range(len(rotatedX)):
        predicted_values.append(ellipse(rotatedX[i],float(results[0]),float(results[1]),float(results[2]),float(results[3]),float(results[4])))
        #print(str(rotatedX[i]) + "\t" + str(predicted_values[i]))
        
    meanSquaredError = np.mean(np.multiply(np.subtract(rotatedY, predicted_values),np.subtract(rotatedY, predicted_values)))
    
    #updates DATA_LIST
    DATA_LIST[INDEX].coefficients.clear() #removes previous entries
    
    for num in results:
        DATA_LIST[INDEX].coefficients.append(num)
    DATA_LIST[INDEX].meanSquaredError = meanSquaredError
    DATA_LIST[INDEX].theta = theta

    
    return meanSquaredError

        
def totalCurvature(a,b, rotatedX0, rotatedXN):
    #curvature of a parabola at some point P = (m,n) is |2a|/(1+(2ax+b)^2)^(3/2), the total curvature 
    # is the integral of the curvature, |a|*(2*a*x+b)/(a*sqrt(4*(a*x)^2 + 4*a*b*x + b^2 + 1))

    X0Curvature = abs(a)*(2*a*rotatedX0+b)/(a*np.sqrt(4*(a*rotatedX0)**2 + 4*a*b*rotatedX0 + b**2 + 1))
    XNCurvature = abs(a)*(2*a*rotatedXN+b)/(a*np.sqrt(4*(a*rotatedXN)**2 + 4*a*b*rotatedXN + b**2 + 1))
    return  abs(float(X0Curvature - XNCurvature)) #total curvature from X0 to XN - the definite integral of curvature


def writeDataToCSV(dataList, outputDirectory):

    outputData = []

    #for each file read, adds the filename and total curvature as a dictionary to outputData
    for i in range(len(dataList)):
        rowData = {"Filename" : dataList[i].filename}

        rotatedX0 = float(rotate(ORIGIN, (dataList[i].X[0], dataList[i].Y[0]), float(dataList[i].theta))[0])
        rotatedXN = float(rotate(ORIGIN, (dataList[i].X[-1], dataList[i].Y[-1]),float(dataList[i].theta))[0])

        rowData["Total Curvature"] = totalCurvature(dataList[i].coefficients[0], dataList[i].coefficients[1], rotatedX0, rotatedXN)
        outputData.append(copy.deepcopy(rowData))
    
    #writes data to csv
    with open(outputDirectory + "TotalCurvatures.csv", "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["Filename","Total Curvature"])
        writer.writeheader()
        writer.writerows(outputData)
        csvfile.close()

def main():

    DATA_LIST = importData(data_path=DATA_DIRECTORY)


    for i in range(len(DATA_LIST)): #prints information related to the optimization
        scipy.optimize.brute(func=rotateCurve,ranges=[slice(0,360,0.1)],args=DATA_LIST[i], Ns=10)
        DATA_LIST[i].print()


        for z in range(len(DATA_LIST[i].X)): #prints the rotated XY coords
            rotatedCoord = rotate(ORIGIN, (DATA_LIST[i].X[z], DATA_LIST[i].Y[z]), DATA_LIST[i].theta)
            print(str(rotatedCoord[0]) + "\t" + str(rotatedCoord[1]))

    writeDataToCSV(DATA_LIST,DATA_DIRECTORY)
    print("Curvatures stored in the .csv in the data Directory")
    
    return 0



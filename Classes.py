import numpy as np
from Imports import *
from Settings import *


#stores information associated with parabolas and parabolic fitting
class Parabola:
    coefficients: dict
    meanAbsolutePercentageError: float
    rotation: float
    rotatedXY: list
    formula = "a*x**2 + b*x + c"
    curvatureFormula = "(2*a*x + b)*Abs(a)/(a*sqrt(4*a**2 * x**2 + 4*a*b*x + b**2 + 1))" # curvature of a parabola in standard form
                                                                                         # at some point P = (m,n) is |2a|/(1+(2ax+b)^2)^(3/2)
    def __init__(self):
        self.coefficients = {}
        self.meanAbsolutePercentageError = -1.0
        self.rotation = 0.0
        self.rotatedXY = []

    #For a given rotation of XY data points, finds the coefficients of the parabola that fits best
    def findRotatedParabola(self, theta:float, filedata) -> float:
        theta = math.radians(theta)
        predicted_values = []

        for i in range(len(filedata.XY)):  # rotates the xy data by a given theta, in radians, probably
            self.rotatedXY.append(rotate(ORIGIN, filedata.XY[i], theta))

        xData = [XY[0] for XY in self.rotatedXY]
        yData = [XY[1] for XY in self.rotatedXY]

        results = scipy.optimize.curve_fit(f=self.standardParabola, xdata=xData, ydata=yData, p0=[1, 1, 1]) #p0 = a, b, c from ax^2+bx+c
        results = results[0] #selects just coefficients

        for x in xData:
            predicted_values.append(self.standardParabola(x, results[0], results[1], results[2]))

        # calculates the mean squared error
        self.meanAbsolutePercentageError = mean_absolute_percentage_error(y_true=yData, y_pred=predicted_values) * 100 #converts to a percentage

        # updates coefficients
        self.coefficients = {'a':results[0], 'b':results[1], 'c':results[2]}

        return self.meanAbsolutePercentageError

    def standardParabola(self, x, a, b, c):
        return a * x ** 2 + b * x + c

    def fitCurve(self, XYList: list, filedata):
        #finds the optimal angle of rotation, theta, about the origin that minimizes the mean squared error associated with
        #the best fitting parabola for the given angle
        self.coefficients['theta'] = scipy.optimize.brute(func=self.findRotatedParabola, ranges=[slice(0, 20, 0.1)], args=filedata, Ns=10)

    def printFormula(self):
        print("\n" + sp.latex(sp.sympify(self.formula).subs({'a':self.coefficients['a'], 'b':self.coefficients['b'], 'c':self.coefficients['c']})))

    def calculateCurvature(self):
        # curvature of a parabola at some point P = (m,n) is |2a|/(1+(2ax+b)^2)^(3/2), the total phragmoplast curvature
        # is the integral of the curvature, |a|*(2*a*x+b)/(a*sqrt(4*(a*x)^2 + 4*a*b*x + b^2 + 1)), over the region the
        # phragmoplast occupies

        a,b,x = sp.symbols("a,b,x")

        X0Curvature = sp.sympify(self.curvatureFormula).subs({a:self.coefficients['a'], b:self.coefficients['b'], x:self.rotatedXY[0][0]})
        XnCurvature = sp.sympify(self.curvatureFormula).subs({a:self.coefficients['a'], b:self.coefficients['b'], x:self.rotatedXY[-1][0]})

        return abs((X0Curvature - XnCurvature).evalf())  # total curvature from X0 to XN - the definite integral of curvature

#stores information associated with ellipses and elliptic fitting
class Ellipse:
    coefficients: dict
    meanAbsolutePercentageError: float
    model = None
    curvatureFormula = "Abs(a*b)/(a**2*cos(t)**2 + b**2*sin(t)**2)**(3/2)"
    formula = ["xc + a*cos(theta)*cos(t) - b*sin(theta)*sin(t)",
               "yc + a*sin(theta)*cos(t) + b*cos(theta)*sin(t)"] #from scikit-image's ellipse model


    def __init__(self):
        self.coefficients = dict()
        self.meanAbsolutePercentageError = -1.0

    def fitCurve(self, XYList:list, filedata): #XYList should be a list of (x,y) tuples
        t = sp.symbols("t")
        XYList = [(np.float64(XY[0]), np.float64(XY[1])) for XY in XYList] #ensures all tuples contain numpy float64's
        XYList = np.array(XYList) #makes XYList usable with scikit-image's EllipseModel



        self.model = skim.EllipseModel
        self.model.estimate(self.model, XYList) #fits the model
        xc, yc, a, b, theta = self.model.params #gets the parameters from the EllipseModel that minimizes mean squared error
                                                #Note: theta represents the angle the ellipse is rotated about it's center,
                                                #as opposed to the meaning of theta used for fitting a parabola

        self.coefficients = {'xc': xc, 'yc': yc, 'a': a, 'b': b, 'theta': theta}

        yTrue = [XY[0] for XY in XYList]
        yPred = []

        X = sp.sympify(self.formula[0]).subs(self.coefficients)
        for XY in XYList:
            T = self.findTClosestToPointOnEllipse(XY)
            yPred.append(X.subs({t:T}).evalf())

        self.meanAbsolutePercentageError = mean_absolute_percentage_error(y_true = yTrue, y_pred = yPred) * 100

    def findCurvatureFormula(self):
        if self.formula != ["xc + a*cos(theta)*cos(t) - b*sin(theta)*sin(t)","yc + a*sin(theta)*cos(t) + b*cos(theta)*sin(t)"]:
            # Finds the general function for curvature for an ellipse of this particular form and prints to the console in latex
            t = sp.symbols("t")
            dX = sp.diff(sp.sympify(self.formula[0]), t)
            dY = sp.diff(sp.sympify(self.formula[1]), t)

            ddX = sp.diff(dX, t)
            ddY = sp.diff(dY, t)

            self.curvatureFormula = (dX * ddY - dY * ddX) / (dX ** 2 + dY ** 2) ** (sp.S(3) / 2)
            self.curvatureFormula = sp.sstr(sp.trigsimp(sp.factor(sp.S(self.curvatureFormula))))  # does some simplification to try to speed up computation

    def curvatureFunction(self, t):
        output = sp.sympify(self.curvatureFormula).subs(self.coefficients).subs('t',t)
        return output.evalf()


    def findTClosestToPointOnEllipse(self, pointOfInterest: tuple) -> float:
        t = sp.symbols("t")
        
        # Set up the parameterized equations of the ellipse being used
        X = sp.sympify(self.formula[0])
        Y = sp.sympify(self.formula[1])

        # Substitute the known values into the equations
        X = X.subs(self.coefficients)
        Y = Y.subs(self.coefficients)

        # makes X & Y lambda functions, for quick evaluation
        X = sp.lambdify(t, X, modules='numpy')
        Y = sp.lambdify(t, Y, modules='numpy')

        distanceFromPointToEllipse = lambda T, x, y: ((X(T) - x) ** 2 + (Y(T) - y) ** 2)

        x0 = np.array([0], dtype=np.float64)
        output = minimize(fun=distanceFromPointToEllipse, x0=x0, method='nelder-mead',args=(pointOfInterest[0],pointOfInterest[1]))

        return output.x[0]


    def calculateCurvature(self, filedata):
        curvature = 0.0

        xc, yc, a, b, theta, m, n = sp.symbols('xc, yc, a, b, theta, m, n', real=True)  # Declare sympy symbols

        self.findCurvatureFormula()

        #finds the limits of integration m,n for the curvature function
        lowerBound = self.findTClosestToPointOnEllipse(filedata.XY[0])
        upperBound = self.findTClosestToPointOnEllipse(filedata.XY[-1])
        
        #update bounds to ensure 0<=t<=2pi, required for Simpson's Rule
        if lowerBound < 0:
            lowerBound = lowerBound + (sp.Abs(sp.floor(lowerBound/(2*sp.pi))) * 2*sp.pi)
        elif lowerBound > 2*sp.pi:
            lowerBound = lowerBound - ((sp.floor(lowerBound/(2*sp.pi))) * 2*sp.pi)


        if upperBound < 0:
            upperBound = upperBound + (sp.Abs(sp.floor(upperBound/(2*sp.pi))) * 2*sp.pi)
        elif lowerBound > 2 * sp.pi:
            lowerBound = upperBound - (sp.floor(upperBound/(2*sp.pi)) * 2*sp.pi)
        
        #ensures upperBound is greater than lowerBound
        if lowerBound > upperBound:
            lowerBound, upperBound = upperBound, lowerBound
        
        
        # Implementation of Simpson's Rule
        numberOfSubIntervals = 100  # arbitrary number, I really rather not implement this in a way that loops until a particular tolerance is met
        stepLength = sp.sympify(sp.S(upperBound - lowerBound)/numberOfSubIntervals)
        xVals = [x for x in np.linspace(start=lowerBound, stop=upperBound, num=numberOfSubIntervals)]
        yVals = [self.curvatureFunction(x) for x in xVals]

        curvature = simpson(y=yVals, x=xVals)
        return curvature

    def printFormula(self):
        print("\n" + "(" + sp.latex(sp.sympify(self.formula[0]).subs(self.coefficients)) + ", " + sp.latex(sp.sympify(self.formula[1]).subs(self.coefficients)) + ")")


# stores data associated with each xy-coordinate file, such as name, the method that provides the best curve fitting, etc
# along with associated functions for displaying or calculating various attributes
class FileData:
    filename: str
    XY: list
    bestFitType: str
    curvature: float
    parabola: Parabola
    ellipse: Ellipse


    def __init__(self, filename="None"):
        self.filename = filename
        self.XY = []
        self.bestFitType = "Undetermined"
        self.curvature = 0.0
        self.parabola = Parabola()
        self.ellipse = Ellipse()

    def updateBestFitType(self):  # updates bestFitType based on which fitted curve has the least mean squared error
        if self.parabola.meanAbsolutePercentageError > 0 and self.ellipse.meanAbsolutePercentageError > 0:
            if self.parabola.meanAbsolutePercentageError > self.ellipse.meanAbsolutePercentageError:
                self.bestFitType = "Ellipse"
            else:
                self.bestFitType = "Parabola"

    def fitCurve(self, curve: Parabola | Ellipse):
        curve.fitCurve(XYList=self.XY, filedata=self)

        self.updateBestFitType()

    def findCurvature(self):

        self.fitCurve(self.ellipse)
        self.fitCurve(self.parabola)

        # calculates the curvature using the curve that has the least mean absolute percentage error
        if self.bestFitType == "Ellipse":
            self.curvature = self.ellipse.calculateCurvature(self)
        elif self.bestFitType == "Parabola":
            self.curvature = self.parabola.calculateCurvature()

    def print(self):
        print("\nfilename: " + self.filename)
        print("\nModel with least square error: " + self.bestFitType)
        print("\nparabolic coefficients: a=" + str(self.parabolicCoefficients[0]) +
              " b=" + str(self.parabolicCoefficients[1]) + " c=" + str(self.parabolicCoefficients[2]) + " theta=" + str(
            self.theta))
        print("\tparabolic mean squared error: " + str(self.parabolicMeanSquaredError))
        print("\nelliptical coefficients: xc=" + str(self.ellipticalCoefficients[0]) +
              " yc=" + str(self.ellipticalCoefficients[1]) + " a=" + str(self.ellipticalCoefficients[2]) +
              " b=" + str(self.ellipticalCoefficients[3]) + " theta=" + str(self.ellipticalCoefficients[4]))
        print("\telliptical mean squared error: " + str(self.ellipticalMeanSquaredError))

    def printXY(self):
        for i in range(len(self.XY)):
            print(str(self.XY[i][0]) + "\t" + str(self.XY[i][1]))

    def printRotatedCoords(self, angle):
        for i in range(len(self.XY)):
            temp = rotate(ORIGIN, self.XY[i], angle)
            print(str(temp[0]) + "\t" + str(temp[1]))




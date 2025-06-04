from __future__ import annotations

from Imports import *

#stores information associated with Bézier and Bézier curve fitting
class Bezier:
    error_info: dict
    control_points: np.ndarray
    curve: dict
    counter = 0

    #loop until some tolerance is met
        #find the ideal locations of control points

    def __init__(self):
        self.error_info = dict()
        self.curve = dict()
    
    # write if necessary for initial_curve_guess
    # def rationa_bezier_to_barycentric(self):
    #     return 0
    # 
    # def barycentric_to_rational_bezier(self):
    #     return 0
    def normalized_absolute_residue_sum(self, y_true=np.ndarray, y_pred=np.ndarray) -> float:
        error = abs(y_true - y_pred)
        error = np.sum(error)
        error = error / np.sum(y_true)
        return error

    def normalized_residue_sum(self, y_true=np.ndarray, y_pred=np.ndarray) -> float:
        error = y_true - y_pred
        error = np.sum(error)
        error = error / np.sum(y_true)
        return error

    def mean_absolute_percent_error(self, y_true:np.ndarray, y_pred:np.ndarray) -> float:

        error = (y_true - y_pred) / y_true
        error = np.sum(abs(error))
        error /= len(y_true)

        return float(error * 100)

    def mean_absolute_log_error(self, y_true:np.ndarray, y_pred:np.ndarray) -> float:

        error = y_true / y_pred
        error = np.log10(error)
        error = np.sum(abs(error))
        error /= len(y_true)

        return float(error)

    def mean_squared_error(self, y_true:np.ndarray, y_pred:np.ndarray) -> float:
        if(len(y_true) != len(y_pred)):
            raise Exception("inputs are of different sizes")

        error = y_true - y_pred
        error = error**2
        error = np.sum(error)
        error /= len(y_true)

        return float(error)

    #use Barycentric form of a Rational Beziér curve to force the curve to fit to the 'vertex' of the phragmoplast
    def initial_curve_guess(self, xy_data):
        interpolation_points = []
        x_expression = 0
        y_expression = 0
        divisor = 0
        weight = 1
        sig = 1

        t = sp.symbols("t", real=True)  # parameterized value

        #finds the barycentric form of a formula that fits the first, middle, and last xy points in xy data
        x_points = [xy_data[0][0], xy_data[int(len(xy_data)/2)][0], xy_data[-1][0]]
        y_points = [xy_data[0][1], xy_data[int(len(xy_data)/2)][1], xy_data[-1][1]]
        t_values = [0, 0.5, 1]

        for i in range(3):
            x_expression += sig * 1 / (t - t_values[i]) * x_points[i]
            y_expression += sig * 1 / (t - t_values[i]) * y_points[i]
            divisor += sig * 1 / (t - t_values[i])
            interpolation_points += [x_points[i], y_points[i]]
            sig = -sig

        x_expression /= divisor
        y_expression /= divisor

        return interpolation_points, x_expression, y_expression, [1,1,1], t_values

    def barycentric_expression(self, interpolation_points:list[tuple], weights:list[float], t_values:list[float]) -> tuple:
        num_interpolation_points = len(interpolation_points)
        x_expression = 0
        y_expression = 0
        common_term = 0
        divisor = 0
        sign = 1

        j = sp.symbols("j", real=True)  # parameterized value

        for i in range(num_interpolation_points):
            common_term = sign * sp.S(weights[i]) / (sp.S(j) - t_values[i])
            x_expression += common_term * interpolation_points[i][0]
            y_expression += common_term * interpolation_points[i][1]
            divisor += common_term
            sign = -sign

        x_expression /= divisor
        y_expression /= divisor

        return str(x_expression), str(y_expression)

    #generates sympy formulas for evaluating a rational Bézier curve
    def rational_Bezier_expression(self, num_control_points:int, control_points:list[tuple], weights:list[float], do_rational=True) -> tuple:
        """ Returns the rational bezier expression for x and y in terms of t as sympy expressions """

        xExpression = 0
        yExpression = 0
        divisor = 0
        t = sp.symbols("t", real=True) #parameterized value


        #implementation of a rational Bézier curve using De Casteljau's algorithm, a closed form solution for calculating Bézier curves of an arbitrary degree
        current_binomial_term = 0
        for i in range(num_control_points):
            current_binomial_term = binom(num_control_points - 1, i)
            divisor +=  weights[i] * current_binomial_term * (1-t)**(num_control_points - i - 1) * t**i
            xExpression += weights[i] * control_points[i][0] * current_binomial_term * (1-t)**(num_control_points - i - 1) * t**i
            yExpression += weights[i] * control_points[i][1] * current_binomial_term * (1-t)**(num_control_points - i - 1) * t**i
        if(do_rational):
            xExpression, yExpression = (xExpression/divisor), (yExpression/divisor)


        return xExpression, yExpression

    #when called, the last arg must be an instance of FileData
    def curve_error(self, *args):
        default_t = 0.499999
        warnings.filterwarnings("ignore")
        control_points, weights = [], []
        XYList = args[1]["XY"]
        t_values = args[1]["t_values"]
        error_measure_method = args[1]["error_measure_method"]
        t = sp.symbols("t", real=True, positive=True)  # declare sympy symbols in order to use sympy


        #sorts data from args[0] into weights and control points
        num_control_points = int(len(args[0])/3)
        weights = args[0][:num_control_points]
        control_points = list(zip(args[0][num_control_points::2],args[0][num_control_points+1::2]))

        #creates rational bezier curves x(t) and y(t) as sympy lambda functions for quick evaluation when measuring curve error
        x_curve, y_curve = self.barycentric_expression(interpolation_points=control_points, weights=weights, t_values=t_values)
        x_curve = sp.sympify(x_curve).subs('j', t)
        y_curve = sp.sympify(y_curve).subs('j', t)

        #lambdify expressions for quick evaluation
        x_curve_lambda, y_curve_lambda = sp.lambdify(t, x_curve, modules="numpy"), sp.lambdify(t, y_curve, modules="numpy")

        def distance_from_bezier(t:float, xCoord:float, yCoord:float) -> float:
            nonlocal x_curve_lambda, y_curve_lambda

            if isinstance(t, np.ndarray):
                t = t[0]

            output = sp.sqrt((x_curve_lambda(t)-xCoord)**2 + (y_curve_lambda(t)-yCoord)**2).evalf()
            return output


        #appends the true y value and the y value predicted from the Bézier curve to yTrue and yPred, respectively, to calculate error
        #jaxxed_x_curve = sympy2jax.SymbolicModule(x_curve)
        #dx_dt = jax.grad(fun=x_curve_lambda)
        #dy_dt = jax.grad(fun=y_curve_lambda)
        dx_dt = sp.lambdify(t, x_curve.diff(t), modules="numpy")
        dy_dt = sp.lambdify(t, y_curve.diff(t), modules="numpy")

        def dd_dt(t:list[np.float64], *args) -> list:
            nonlocal x_curve_lambda, y_curve_lambda
            out = []
            for i in t:
                out.append(((x_curve_lambda(i)-args[0]) * dx_dt(i) + (y_curve_lambda(i) - args[1]) * dy_dt(i)))
            if type(out) == type(None): out = 0
            return out

        true_y, predicted_y = np.empty(0), np.empty(0)
        for XY in XYList:
            true_y = np.append(true_y, XY[1])
            closest_t = scipy_optimize.fsolve(func=dd_dt, x0=default_t, args=XY)

            if closest_t == default_t:
                closest_t = scipy_optimize.fsolve(func=dd_dt, x0=default_t + 0.1, args=XY)

            predicted_y = np.append(predicted_y, y_curve_lambda(closest_t[0]))

        #selects which error measure method to use based on error_measure_method
        if error_measure_method == "mean_squared_error": error = self.mean_squared_error(y_true=true_y, y_pred=predicted_y)
        elif error_measure_method == "mean_absolute_log_error":error = self.mean_absolute_log_error(y_true=true_y, y_pred=predicted_y)
        elif error_measure_method == "mean_absolute_percent_error": error = self.mean_absolute_percent_error(y_true=true_y, y_pred=predicted_y)
        elif error_measure_method == "normalized_residue_sum": error = self.normalized_residue_sum(y_true=true_y, y_pred=predicted_y)
        elif error_measure_method == "normalized_absolute_residue_sum": error = self.normalized_absolute_residue_sum(y_true=true_y, y_pred=predicted_y)
        else: raise Exception(f"Unknown error measure method: {error_measure_method}")

        self.error_info = {"error":error, "true_y":true_y, "predicted_y":predicted_y}
        return error

    #fits a rational Bézier curve to the data set by optimizing control points, control point weights, and
    # elevating the degree of the curve as necessary
    def fit_curve(self, filedata:FileData):
        tolerance = 0.00001
        error = 100
        iterationCounter = 1
        max_iterations = 3
        options = {"maxiter":10}
        curvature = 0

        control_points, self.curve["x_curve"], self.curve["y_curve"], control_weights, t_values = self.initial_curve_guess(filedata.XY)
        self.curve["x_curve"], self.curve["y_curve"] = str(self.curve["x_curve"]), str(self.curve["y_curve"])
        num_control_points = int(len(control_points)/2)
        args_curve_error = {"XY": filedata.XY, "t_values": t_values, "error_measure_method": "mean_absolute_log_error"}

        old_error = self.curve_error(control_weights + control_points, args_curve_error)

        weight_bounds = (0, 5)
        coord_bound = (1, None)
        bounds = [weight_bounds, weight_bounds, weight_bounds, coord_bound, coord_bound, coord_bound, coord_bound, coord_bound, coord_bound]
        x0 = control_weights + control_points

        while True: #iteratively refines the barycentric form of a rational Bézier curve by adding more control points until error falls below tolerance or max_iterations is met
            #does curve fitting
            control_points = scipy_optimize.minimize(method='Nelder-Mead', fun=self.curve_error, x0=x0, bounds=bounds, options=options, args=args_curve_error).x.tolist()

            #separates the control weights and control points from scipy optimization
            control_weights = control_points[:num_control_points]
            control_points = control_points[num_control_points:]

            #saves the current curve into self.curve, for calculating curvature
            x_curve, y_curve = self.barycentric_expression(interpolation_points=list(zip(control_points[::2],control_points[1::2])), weights=control_weights, t_values=t_values)
            self.curve["x_curve"], self.curve["y_curve"] = x_curve.replace('j','t'), y_curve.replace('j','t')

            error = self.curve_error(control_weights + control_points, args_curve_error) #error of new curve

            #checks if conditions are met to exit loop
            if error < tolerance or iterationCounter >= max_iterations:
                break
            iterationCounter += 1

            #need to elevate curve
            temp = (abs(self.error_info["predicted_y"] - self.error_info["true_y"]))
            new_interpolation_point = self.new_interpolation_point(filedata.XY, temp, t_values)
            control_points, control_weights, t_values = self.elevate_barycentric_curve(control_points, control_weights, t_values, new_interpolation_point)
            num_control_points = int(len(control_points) / 2)

        curvature = self.calculate_curvature(x_curve, y_curve, t_values)
        return curvature

    def new_interpolation_point(self, xy:list[floats], abs_residues:list[floats], t_values:list[float]) -> dict:
        """Finds the x, y, and t for the point on the barycentric curve with the highest error"""
        small_val = 0.0001
        default_t = 0.5 - small_val
        t = sp.symbols("t", real=True)
        new_point = {"x":-1, "y":-1, "t":-1}
        bound = scipy_optimize.Bounds(0)


        index_of_most_error = np.argmax(abs_residues)
        if type(index_of_most_error) != np.int64:  #incase np.argmax returns an ndarray
            index_of_most_error = index_of_most_error[0]

        #get x_curve and y_curve and turn into a sympy expression
        x_curve, y_curve = sp.parse_expr(self.curve["x_curve"],local_dict={'t':t}), sp.parse_expr(self.curve["y_curve"],local_dict={'t':t})

        dx_dt = sp.lambdify(t, x_curve.diff(t), modules="numpy")
        dy_dt = sp.lambdify(t, y_curve.diff(t), modules="numpy")

        x_curve_lambda, y_curve_lambda = sp.lambdify(t, x_curve, modules="numpy"), sp.lambdify(t, y_curve, modules="numpy")

        def dd_dt(k:list[np.float64], *args) -> list:
            nonlocal x_curve_lambda, y_curve_lambda
            out = []
            for i in k:
                out.append(((x_curve_lambda(i) - args[0]) * dx_dt(i) + (y_curve_lambda(i) - args[1]) * dy_dt(i)))
            if type(out) == type(None): out = 0
            return out

        #most likely to break
        is_close = []
        closest_t, i = 0, 0
        while True:
            closest_t = scipy_optimize.minimize(fun=dd_dt, x0=t_values[i]+small_val, args=xy[index_of_most_error], bounds=bound, method="Nelder-Mead").x[0]

            is_close = [math.isclose(closest_t, t, abs_tol = small_val) for t in t_values] #determines if closest_t is around one of the singularities that occur at t_vals
            if not True in is_close: #if closest_t is not near a singularity, it is probably is correct, so break the loop
                break
            i += 1

        new_point["t"] = closest_t
        new_point["x"], new_point["y"] = x_curve_lambda(closest_t), y_curve_lambda(closest_t)
        return new_point

    #Implementation of proposition 7 from Ramanantoanina and Hormann's 2021 paper "New shape control tools for rational Bézier curve design"
    def elevate_barycentric_curve(self, interpolation_points:list[float], weights:list[float], t_values:list[float], new_interpolation_point_values:dict):
        new_interpolation_point_index = -1
        new_weights = []
        new_t_values = copy.deepcopy(t_values)

        #updates interpolation_points and t_values to include new interpolation point
        for i in range(len(new_t_values)):
            if new_interpolation_point_values["t"] < new_t_values[i]:
                new_t_values.insert(i, new_interpolation_point_values["t"])
                interpolation_points.insert(i, (new_interpolation_point_values["x"], new_interpolation_point_values["y"]))
                new_interpolation_point_index = i
                break

        if new_interpolation_point_index == -1:
            raise Exception("Could not find suitable index for new interpolation point")

        #updates weights to compensate for the new interpolation point,
        for n in range(len(new_t_values)):
            if n < new_interpolation_point_index:
                new_weights.append(weights[n] / (new_interpolation_point_values["t"] - t_values[n]))
            elif n > new_interpolation_point_index:
                new_weights.append(weights[n - 1] / (t_values[n - 1] - new_interpolation_point_values["t"]))
            else: #when n == new_interpolation_point_index
                temp_weight = 0
                sign = (-1) ** (len(weights) + new_interpolation_point_index + 1)
                for j in range(len(weights)):
                        temp_weight +=  sign * weights[j] / (t_values[j] - new_interpolation_point_values["t"])
                        sign = -sign
                new_weights.append(temp_weight)

        return interpolation_points, new_weights, new_t_values

    def calculate_curvature(self, x_curve:str, y_curve:str, t_values:list[np.float64]) -> tuple:
        """returns the sum of curvature over a parametric curve from 0 to 1, returns the calculated curvature and the error associated"""
        t = sp.symbols("t", real=True)
        small_val = sp.S(0.00000000001)
        curvature = 0

        for i in range(len(t_values)):
            t_values[i] = np.float64(t_values[i])

        #ensures x_curve & y_curve are copied appropriately
        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('j','t'),local_dict={'t':t}), sp.parse_expr(str(y_curve).replace('j','t'),local_dict={'t':t})

        dx_dt, dy_dt = x_curve.diff(t), y_curve.diff(t) #finding the 1st derivatives
        ddx_dt, ddy_dt = sp.diff(dx_dt, t), sp.diff(dy_dt, t) #finding the 2nd derivatives

        numerator = abs(dx_dt * ddy_dt - dy_dt * ddx_dt)
        denominator = (dx_dt * dx_dt + dy_dt * dy_dt)**sp.S(3/2)
        integrand = sp.lambdify(t, numerator / denominator,modules="numpy")
        curvature = quad(integrand, 0, 1, points=t_values) #approximates the integral with little error, much faster than using sympy

        return curvature[0], curvature[1]
# stores data associated with each xy-coordinate file, such as name, the method that provides the best curve fitting, etc
# along with associated functions for displaying or calculating various attributes
class FileData:
    filename: str
    XY: list
    curvature: float
    bezier: Bezier


    def __init__(self, filename="None"):
        self.filename = filename
        self.XY = []
        #self.bestFitType = "Undetermined"
        self.curvature = 0.0
        self.bezier = Bezier()
        #self.parabola = Parabola()
        #self.ellipse = Ellipse()

    def findCurvature(self) -> float:
        self.curvature = self.bezier.fit_curve(self)
        return self.curvature

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

class Tests:
    def All_Tests(self):

        print("Barycentric elevation test: " + str(self.Barycentric_Elevation_Test()))
        print("Curvature calculation test: " + str(self.Curvature_Test()))
        return 0

    def Barycentric_Elevation_Test(self) -> bool:
        # testing degree elevation
        test = Bezier()
        output = False
        interpolation_points = [(1, 0), (12 / 13, 5 / 13), (0, 1)]
        t_vals = [0, 1 / 3, 1]
        weights = [2, 13 / 6, 1 / 2]

        x_curve, y_curve = test.barycentric_expression(interpolation_points, weights, t_vals)
        print(sp.latex(sp.S(x_curve.replace('j', 't'))) + "\n" + sp.latex(sp.S(y_curve.replace('j', 't'))))
        new_point = {"x": 3 / 5, "y": 4 / 5, "t": 2 / 3}

        expected_interpolation_points = [(1.0, 0.0), (12.0 / 13.0, 5.0 / 13.0), (3.0 / 5.0, 4.0 / 5.0), (0.0, 1.0)]
        expected_t_vals = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
        expected_weights = [3.0, 13.0 / 2.0, 5.0, 1.4999999999999998]

        interpolation_points, weights, t_vals = test.elevate_barycentric_curve(interpolation_points, weights, t_vals, new_point)
        interpolation_points, weights, t_vals = [(float(x[0]),float(x[1])) for x in interpolation_points], [float(x) for x in weights], [float(x) for x in t_vals]

        if expected_interpolation_points == interpolation_points and expected_t_vals == t_vals and expected_weights == weights:
            output = True

        return output

    def Curvature_Test(self):
        test = Bezier()
        output = False
        interpolation_points = [(1, 0), (sp.S(12) / 13, sp.S(5) / 13), (0, 1)]
        t_vals = [0, sp.S(1) / 3, 1]
        weights = [2, sp.S(13) / 6, sp.S(1) / 2]

        x_curve, y_curve = test.barycentric_expression(interpolation_points, weights, t_vals)
        curvature = test.calculate_curvature(x_curve, y_curve, t_vals)[0]

        if math.isclose(curvature, 1, abs_tol=0.001):
            output = True


        return output
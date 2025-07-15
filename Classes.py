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

    def mean_squared_error(self, y_true:np.ndarray, y_pred:np.ndarray) -> float:
        if(len(y_true) != len(y_pred)):
            raise Exception("inputs are of different sizes")

        error = y_true - y_pred
        error = error**2
        error = np.sum(error)
        error /= len(y_true)

        return float(error)

    def sum_squared_distance(self, x_true:np.ndarray, y_true:np.ndarray, x_pred:np.ndarray, y_pred:np.ndarray) -> float:
        error = 0

        for i in range(len(x_true)):
            error += sp.sqrt((x_true[i] - x_pred[i]) ** 2 + (y_true[i] - y_pred[i]) ** 2)

        return error

    def rotate_point(self, origin, point, angle) -> tuple:
        # Code adapted to use sympy from Mark Dickinson's post, found at
        # https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python

        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """

        ox, oy = origin
        px, py = point

        cos, sin = sp.cos(angle), sp.sin(angle)
        px_ox, py_oy = px - ox, py - oy
        qx = ox + cos * px_ox  - sin * py_oy
        qy = oy + sin * px_ox + cos * py_oy

        return qx, qy

    def rotate_point_set(self, origin, points:list[tuple], angle) -> list[tuple]:
        # Code adapted to use sympy and to accept a list of points from Mark Dickinson's post, found at
        # https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python

        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """

        output = []
        ox, oy = origin
        for point in points:

            px, py = point
            cos, sin = sp.cos(angle), sp.sin(angle)
            px_ox, py_oy = px - ox, py - oy
            qx = sp.S(ox + cos * px_ox - sin * py_oy)
            qy = sp.S(oy + sin * px_ox + cos * py_oy)

            output.append((qx, qy))

        return output

    def scale_point_set(self, points:list[tuple], scaling_factor:float) -> list[tuple]:
        output = []
        for point in points: output.append((point[0] * scaling_factor, point[1] * scaling_factor))
        return output

    def translate_point_set(self, points:list[tuple], horizontal_shift:float, vertical_shift:float) -> list[tuple]:
        """Moves points in list to (x-horizontal_shift, y+vertical_shift)"""

        output = []
        for point in points:
            output.append((point[0] + horizontal_shift, point[1] + vertical_shift))
        return output

    #rotates, scales, and translates data to be concave down and on the provided domain
    def normalize_data(self, xy_data:list[tuple], domain:tuple[float, float]) -> tuple:

        #rotates data so that the p0=(x0, y0) & pN=(x0 , ?)
        center_of_rotation = (xy_data[0][0], xy_data[0][1])
        is_oriented_top_to_bottom = xy_data[1][1] - xy_data[0][1] < 0

        # if math.isclose(xy_data[0][0], xy_data[-1][0]):
        #     theta = sp.pi / 2
        # elif math.isclose(xy_data[0][1], xy_data[-1][1]):
        #     theta = 0
        # elif is_oriented_top_to_bottom:
        #     theta = sp.pi / 2 - sp.atan((xy_data[-1][1] - xy_data[0][1]) / (sp.S(xy_data[-1][0]) - xy_data[0][0]))
        # else:
        #     theta = -1 * sp.atan((xy_data[-1][1] - xy_data[0][1]) / (sp.S(xy_data[-1][0]) - xy_data[0][0]))
        theta = -1 * sp.atan((xy_data[-1][1] - xy_data[0][1]) / (sp.S(xy_data[-1][0]) - xy_data[0][0]))
        new_xy_data = self.rotate_point_set(center_of_rotation, xy_data, theta)

        is_oriented_vertically = math.isclose(new_xy_data[0][0], new_xy_data[-1][0])

        #ensures data is either concave up or down
        if is_oriented_vertically:
            new_xy_data = self.rotate_point_set(center_of_rotation, xy_data, sp.pi/2)
            theta += sp.pi/2

        #finds the index of the middlemost point
        if new_xy_data[0][0] > new_xy_data[-1][0]:
            midpoint_x = new_xy_data[0][0] - (new_xy_data[0][0] -  new_xy_data[-1][0])/2
        else:
            midpoint_x = new_xy_data[-1][1] - (new_xy_data[-1][1] - new_xy_data[0][0])/2

        #finds the index of the tuple with an x-value closest to the midpoint
        closest_index_to_midpoint = np.abs(np.array([i[0] for i in new_xy_data[1:]]) - midpoint_x).argmin() + 1

        #finds if the data is oriented concave down
        is_concave_down = False
        is_oriented_vertically = math.isclose(new_xy_data[0][0], new_xy_data[-1][0])
        is_left_to_right = new_xy_data[0][0] < new_xy_data[-1][0]

        if not is_oriented_vertically:
            if is_left_to_right:
                is_positive_slope = (new_xy_data[closest_index_to_midpoint][1] - new_xy_data[0][1]) / (new_xy_data[closest_index_to_midpoint][0] - new_xy_data[0][0]) >= 0
            else:
                is_positive_slope = (new_xy_data[closest_index_to_midpoint][1] - new_xy_data[0][1]) / (new_xy_data[closest_index_to_midpoint][0] - new_xy_data[0][0]) <= 0

            if is_positive_slope: is_concave_down = True

        #rotates data to be concave down if necessary
        if not is_concave_down:
            new_xy_data = self.rotate_point_set(center_of_rotation, new_xy_data, sp.pi)
            theta += sp.pi

        #scales data so that |xN-x0| is equal to the difference of x-values in the given domain
        scale_factor = abs((domain[1] - domain[0]) / sp.S(new_xy_data[0][0] - new_xy_data[-1][0]))
        new_xy_data = self.scale_point_set(new_xy_data, scale_factor)

        #scales data so that the difference between most extreme x values fits in the domain
        new_x_vals_temp = np.array([i[0] for i in new_xy_data])
        min_x_index, max_x_index = new_x_vals_temp.argmin(), new_x_vals_temp.argmax()

        if not math.isclose(new_xy_data[max_x_index][0] - new_xy_data[min_x_index][0], domain[1] - domain[0]):
            additional_scale_factor = abs((domain[1] - domain[0]) / sp.S(new_xy_data[max_x_index][0] - new_xy_data[min_x_index][0]))
            new_xy_data = self.scale_point_set(new_xy_data, additional_scale_factor)
            scale_factor = scale_factor * additional_scale_factor

        #translates data so that xy-values are on the provided domain
        horizontal_shift, vertical_shift = domain[0] - new_xy_data[min_x_index][0], 0 - new_xy_data[0][1]
        new_xy_data = self.translate_point_set(new_xy_data, horizontal_shift, vertical_shift)

        # for i in range(len(xy_data)):
        #     print(str(new_xy_data[i][0]) + "\t" + str(new_xy_data[i][1]))

        return new_xy_data, {"points":new_xy_data, "theta":theta, "scale_factor":scale_factor, "horizontal_shift":horizontal_shift, "vertical_shift":vertical_shift}

        # rotates, scales, and translates data to be concave down and on the provided domain

    def normalize_data_without_scaling(self, xy_data: list[tuple], domain: tuple[float, float]) -> tuple:

        # rotates data so that the p0=(x0, y0) & pN=(x0 , ?)
        center_of_rotation = (xy_data[0][0], xy_data[0][1])
        is_oriented_top_to_bottom = xy_data[1][1] - xy_data[0][1] < 0

        if math.isclose(xy_data[0][0], xy_data[-1][0]):
            theta = sp.pi / 2
        elif math.isclose(xy_data[0][1], xy_data[-1][1]):
            theta = 0
        elif is_oriented_top_to_bottom:
            theta = sp.pi / 2 - sp.atan((xy_data[-1][1] - xy_data[0][1]) / (sp.S(xy_data[-1][0]) - xy_data[0][0]))
        else:
            theta = -1 * sp.atan((xy_data[-1][1] - xy_data[0][1]) / (sp.S(xy_data[-1][0]) - xy_data[0][0]))
        new_xy_data = self.rotate_point_set(center_of_rotation, xy_data, theta)

        is_oriented_vertically = math.isclose(new_xy_data[0][0], new_xy_data[-1][0])

        # finds the index in the data set that's x or y (depending on if it is oriented vertically) lies closest to the middle x or y
        if is_oriented_vertically:
            if new_xy_data[0][1] > new_xy_data[-1][1]:
                midpoint_y = new_xy_data[0][1] - (new_xy_data[0][1] - new_xy_data[-1][1]) / 2
            else:
                midpoint_y = new_xy_data[-1][1] - (new_xy_data[-1][1] - new_xy_data[0][1]) / 2
            # finds the index of the tuple with a y-value closest to the midpoint
            closest_index_to_midpoint = np.abs(np.array([i[1] for i in new_xy_data[1:]]) - midpoint_y).argmin() + 1
        else:
            if new_xy_data[0][0] > new_xy_data[-1][0]:
                midpoint_x = new_xy_data[0][0] - (new_xy_data[0][0] - new_xy_data[-1][0]) / 2
            else:
                midpoint_x = new_xy_data[-1][1] - (new_xy_data[-1][1] - new_xy_data[0][0]) / 2
            # finds the index of the tuple with an x-value closest to the midpoint
            closest_index_to_midpoint = np.abs(np.array([i[0] for i in new_xy_data[1:]]) - midpoint_x).argmin() + 1

        # update/define values for rotating data points to be concave down, if necessary

        is_positive_slope = (new_xy_data[closest_index_to_midpoint][1] - new_xy_data[0][1]) / (
                    new_xy_data[closest_index_to_midpoint][0] - new_xy_data[0][0]) > 0
        is_oriented_top_to_bottom = new_xy_data[1][1] - new_xy_data[0][1] < 0

        # check if data is concave down
        if is_oriented_vertically == False and new_xy_data[int(len(new_xy_data) / 2)][1] > new_xy_data[0][1]:
            is_concave_down = True
        else:
            is_concave_down = False

        # rotates data to be concave down data is not concave down
        if not is_concave_down:
            if is_positive_slope == False and is_oriented_top_to_bottom == False:
                new_xy_data = self.rotate_point_set(center_of_rotation, new_xy_data, -sp.pi / 2)
                theta = theta - sp.pi / 2
            elif is_positive_slope == False and is_oriented_top_to_bottom:
                new_xy_data = self.rotate_point_set(center_of_rotation, new_xy_data, sp.pi / 2)
                theta += sp.pi / 2
            elif is_positive_slope and is_oriented_top_to_bottom == False:
                new_xy_data = self.rotate_point_set(center_of_rotation, new_xy_data, sp.pi / 2)
                theta += sp.pi / 2
            elif is_positive_slope and is_oriented_top_to_bottom:
                new_xy_data = self.rotate_point_set(center_of_rotation, new_xy_data, -sp.pi / 2)
                theta = theta - sp.pi / 2

        # # scales data so that |xN-x0| is equal to the difference of x-values in the given domain
        # scale_factor = abs((domain[1] - domain[0]) / sp.S(new_xy_data[0][0] - new_xy_data[-1][0]))
        # new_xy_data = self.scale_point_set(new_xy_data, scale_factor)

        #finds most extreme x_vals
        new_x_vals_temp = np.array([i[0] for i in new_xy_data])
        min_x_index, max_x_index = new_x_vals_temp.argmin(), new_x_vals_temp.argmax()

        ## scales data so that the difference between most extreme x values fits in the domain
        # if not math.isclose(new_xy_data[max_x_index][0] - new_xy_data[min_x_index][0], domain[1] - domain[0]):
        #     additional_scale_factor = abs(
        #         (domain[1] - domain[0]) / sp.S(new_xy_data[max_x_index][0] - new_xy_data[min_x_index][0]))
        #     new_xy_data = self.scale_point_set(new_xy_data, additional_scale_factor)
        #     scale_factor = scale_factor * additional_scale_factor

        # translates data so that xy-values are on the provided domain
        horizontal_shift, vertical_shift = domain[0] - new_xy_data[min_x_index][0], 0 - new_xy_data[0][1]
        new_xy_data = self.translate_point_set(new_xy_data, horizontal_shift, vertical_shift)

        # for i in range(len(xy_data)):
        #     print(str(new_xy_data[i][0]) + "\t" + str(new_xy_data[i][1]))

        return new_xy_data, {"points": new_xy_data, "theta": theta,
                             "horizontal_shift": horizontal_shift, "vertical_shift": vertical_shift}

    def denormalize_data(self,points:list[tuple],theta:float,scale_factor:float,horizontal_shift:float,vertical_shift:float) -> list[tuple]:


        new_xy_data = self.translate_point_set(points, -1*horizontal_shift, -1*vertical_shift)
        new_xy_data = self.scale_point_set(new_xy_data, 1/scale_factor)
        center_of_rotation = (new_xy_data[0][0], new_xy_data[0][1])
        new_xy_data = self.rotate_point_set(center_of_rotation, new_xy_data, 2*sp.pi - theta)

        return new_xy_data

    def rotate_curve(self,x_curve:str, y_curve:str, origin:tuple[float,float], t_origin:float, theta:float) -> tuple:
        new_x_curve, new_y_curve  = "", ""
        t = sp.symbols("t", real=True)

        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('t', 't'), local_dict={'t': t}), sp.parse_expr(
            str(y_curve).replace('t', 't'), local_dict={'t': t})
        new_x_curve = origin[0] + sp.cos(theta) * (x_curve - origin[0]) - sp.sin(theta) * (y_curve - origin[1])
        new_y_curve = origin[1] + sp.sin(theta) * (x_curve - origin[0]) + sp.cos(theta) * (y_curve - origin[1])

        return str(new_x_curve), str(new_y_curve)

    def scale_curve(self, x_curve:str, y_curve:str, scale_factor:float):
        new_x_curve, new_y_curve = "", ""
        t = sp.symbols("t", real=True)

        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('t', 't'), local_dict={'t': t}), sp.parse_expr(
            str(y_curve).replace('t', 't'), local_dict={'t': t})

        new_x_curve, new_y_curve = x_curve * scale_factor, y_curve * scale_factor

        return str(new_x_curve), str(new_y_curve)

    def translate_curve(self, x_curve:str, y_curve:str, horizontal_shift:float, vertical_shift:float) -> tuple:
        new_x_curve, new_y_curve = "", ""
        t = sp.symbols("t", real=True)

        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('t', 't'), local_dict={'t': t}), sp.parse_expr(
            str(y_curve).replace('t', 't'), local_dict={'t': t})

        new_x_curve, new_y_curve = x_curve + horizontal_shift, y_curve + vertical_shift
        return str(new_x_curve), str(new_y_curve)

    #should only work when (x(t0), y(t0)) ~= p0
    def denormalize_curve(self, points:list[tuple], origin:tuple[float,float], x_curve:str, y_curve:str, theta:float, scale_factor:float, horizontal_shift:float, vertical_shift:float) -> tuple:
        small_val = 0.00000001
        new_x_curve, new_y_curve = "", ""
        t = sp.symbols("t", real=True)

        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('j', 't'), local_dict={'t': t}), sp.parse_expr(
            str(y_curve).replace('j', 't'), local_dict={'t': t})

        new_x_curve, new_y_curve = self.translate_curve(x_curve, y_curve, -horizontal_shift, -vertical_shift)
        new_x_curve, new_y_curve = self.scale_curve(new_x_curve, new_y_curve, 1 / scale_factor)

        def dist_due_to_theta(t: list[np.float64], *args) -> list:
            nonlocal points, origin
            out, temp, dist =  [], [], []
            for i in t:
                dist = []
                temp = self.rotate_point_set(origin,points,i)
                for i in range(len(temp)):
                    dist.append(sqeuclidean(temp[i], points[i]))
                dist = sum(dist)
                out.append(dist)

            return out

        theta = scipy_optimize.fsolve(func=dist_due_to_theta, x0=theta, args=[])[0]

        rot_x_curve, rot_y_curve = self.rotate_curve(new_x_curve, new_y_curve, origin, small_val, -theta)


        return str(new_x_curve), str(new_y_curve)


    def initial_curve_guess(self, xy_data:list[tuple], curve_error_method:str="sum_squared_distance") -> tuple:
        interpolation_points, interpolation_weights, old_weights, errors = [], [], [], []
        t_values = [0, 0.5, 1]
        args_curve_error = {"XY":xy_data, "t_values": t_values, "error_measure_method": curve_error_method}
        small_val = 0.000001
        options = {"maxiter": 5}

        #finds the distance between the first point and every other non-endpoint and
        # the distance between the last point and every other non-endpoint
        distance = [[],[]]
        for i in range(1, len(xy_data) - 1):
            distance[0].append(sqeuclidean(xy_data[i], xy_data[-1]))
            distance[1].append(sqeuclidean(xy_data[i], xy_data[0]))

        #converts distance to store the absolute differences between the distances,
        # used to find the middle-ish most point, as it is not guaranteed that the points are spaced equally
        distance = abs(np.array(distance[0]) - np.array(distance[1]))

        middle_most_point = np.argmin(distance)
        # x_vals_only = np.array([i[0] for i in xy_data])
        # min_x_index, max_x_index = x_vals_only.argmin(), x_vals_only.argmax()
        interpolation_points = [xy_data[0][0], xy_data[0][1], xy_data[middle_most_point][0],
                                xy_data[middle_most_point][1], xy_data[-1][0], xy_data[-1][1]]
        interpolation_weights = [1, 5, 1]

        return interpolation_points, interpolation_weights, t_values

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
        small_val = 0.000001
        default_t = 0.5 - small_val
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

        #appends the true y value and the y value predicted from the Bézier curve to yTrue and yPred, respectively, to calculate error
        dx_dt = sp.lambdify(t, x_curve.diff(t), modules="numpy")
        dy_dt = sp.lambdify(t, y_curve.diff(t), modules="numpy")

        #derivative of ~distance from XY with respect to t
        def dd_dt(t:list[np.float64], *args) -> list:
            nonlocal x_curve_lambda, y_curve_lambda, dx_dt, dy_dt
            out = []
            for i in t:
                out.append(((x_curve_lambda(i)-args[0]) * dx_dt(i) + (y_curve_lambda(i) - args[1]) * dy_dt(i)))
            return out

        #finds the value for t closest to each point, then finds that point on the provided curve
        closest_t = 0
        true_x, predicted_x, true_y, predicted_y = np.empty(0), np.empty(0), np.empty(0), np.empty(0)

        if error_measure_method == "mean_squared_error":

            def objective_function(t:list[np.float64], *args) -> list:
                nonlocal x_curve_lambda
                out = []
                for i in t:
                    out.append((x_curve_lambda(i)-args[0])**2)
                return out

            for i in range(len(XYList)):
                true_x, true_y = np.append(true_x, XYList[i][0]), np.append(true_y, XYList[i][1])
                t1 = scipy_optimize.fsolve(func=objective_function, x0=small_val, args=XYList[i])[0]
                t2 = scipy_optimize.fsolve(func=objective_function, x0=1 - small_val, args=XYList[i])[0]

                # ensures 0<=t1,t2<=1
                t1, t2 = min(max(0 + small_val, t1), 1 - small_val), min(max(0 + small_val, t2), 1 - small_val)

                if (x_curve_lambda(t1) -  XYList[i][0])**2 < (x_curve_lambda(t2) -  XYList[i][0])**2:
                    closest_t = t1
                else:
                    closest_t = t2

                predicted_x, predicted_y = np.append(predicted_x, x_curve_lambda(closest_t)), np.append(predicted_y,y_curve_lambda(closest_t))


        else:
            for i in range(len(XYList)):
                true_x, true_y = np.append(true_x, XYList[i][0]), np.append(true_y, XYList[i][1])
                t1 = scipy_optimize.fsolve(func=dd_dt, x0=small_val, args=XYList[i])[0]
                t2 = scipy_optimize.fsolve(func=dd_dt, x0=1 - small_val, args=XYList[i])[0]

                #ensures 0<=t1,t2<=1
                t1, t2 = min(max(0 + small_val, t1), 1 - small_val), min(max(0 + small_val, t2), 1 - small_val)

                if sqeuclidean((x_curve_lambda(t1),y_curve_lambda(t1)), XYList[i]) < sqeuclidean((x_curve_lambda(t2),y_curve_lambda(t2)), XYList[i]):
                    closest_t = t1
                else:
                    closest_t = t2


                predicted_x, predicted_y = np.append(predicted_x, x_curve_lambda(closest_t)), np.append(predicted_y, y_curve_lambda(closest_t))

        #selects which error measure method to use based on error_measure_method 
        if error_measure_method == "mean_squared_error": error = self.mean_squared_error(y_true=true_y, y_pred=predicted_y)
        elif error_measure_method == "sum_squared_distance": error = self.sum_squared_distance(x_true=true_x, y_true=true_y, x_pred=predicted_x, y_pred=predicted_y)
        else: raise Exception(f"Unknown error measure method: {error_measure_method}")

        self.error_info = {"error":copy.deepcopy(error), "true_y":copy.deepcopy(true_y), "predicted_y":copy.deepcopy(predicted_y),"true_x":copy.deepcopy(true_x),"predicted_x":copy.deepcopy(predicted_x),}
        return error

    #fits a rational Bézier curve to the data set by optimizing control points, control point weights, and
    # elevating the degree of the curve as necessary
    def fit_curve(self, filedata:FileData):
        tolerance = 0.00001
        error = 100
        iteration_counter = 1
        max_iterations = 5
        options = {"maxiter":50}
        normalization_domain = (10, 110)
        multiplier = 1.025
        curvature = 0
        t = sp.symbols("t", real=True)

        #applies an affine transformation on the data so that it is (mostly) on the normalization_domain
        # normalized_xy, normalization_info = self.normalize_data(filedata.XY, normalization_domain)
        normalized_xy, normalization_info = self.normalize_data(filedata.XY, normalization_domain)
        for i in normalized_xy: print(str(i[0]) + "\t" + str(i[1]))

        control_points, control_weights, t_values = self.initial_curve_guess(normalized_xy, "mean_squared_error")

        guess_x, guess_y = self.barycentric_expression(list(zip(control_points[::2],control_points[1::2])), control_weights, t_values)
        guess_x, guess_y = sp.sympify(guess_x).subs('j', t), sp.sympify(guess_y).subs('j', t)
        print("Initial guess: " + "(" + sp.latex(sp.S(guess_x)) + "," + sp.latex(sp.S(guess_y)) + ")")

        args_curve_error = {"XY": normalized_xy, "t_values": t_values, "error_measure_method": "mean_squared_error"}

        old_error = self.curve_error(control_weights + control_points, args_curve_error)

        num_control_points = int(len(control_points) / 2)

        x_only, y_only = [i[0] for i in normalized_xy], [i[1] for i in normalized_xy]


        weight_bounds = (0, None)
        x_bound = (min(x_only) / multiplier, max(x_only) * multiplier)
        y_bound = (min(y_only) / multiplier, max(y_only) * multiplier)
        #coord_bound = (1, None)
        bounds = [weight_bounds, weight_bounds, weight_bounds, (control_points[0], control_points[0]), #fixes the endpoints in place
            (control_points[1], control_points[1]), x_bound, y_bound, (control_points[-2], control_points[-2]),
            (control_points[-1], control_points[-1])]


        while True: #iteratively refines the barycentric form of a rational Bézier curve by adding more control points until error falls below tolerance or max_iterations is met
            x0 = control_weights + control_points

            #does curve fitting
            control_points = scipy_optimize.minimize(method='Nelder-Mead', fun=self.curve_error, x0=x0, bounds=bounds, options=options, args=args_curve_error).x.tolist()

            #separates the control weights and control points from scipy optimization
            control_weights = control_points[:num_control_points]
            control_points = control_points[num_control_points:]


            #saves the current curve into self.curve, for calculating curvature
            x_curve, y_curve = self.barycentric_expression(interpolation_points=list(zip(control_points[::2],control_points[1::2])), weights=control_weights, t_values=t_values)
            self.curve["x_curve"], self.curve["y_curve"] = x_curve.replace('j','t'), y_curve.replace('j','t')


            error = self.curve_error(control_weights + control_points, args_curve_error) #error of new curve
            #print(filedata.filename + ": " + str(iteration_counter))
            print("Error: " + str(error) + "\t" + "Curve formula: " + "(" + sp.latex(sp.S(self.curve["x_curve"])) + "," + sp.latex(sp.S(self.curve["y_curve"])) + ")")

            #checks if conditions are met to exit loop
            if error < tolerance or iteration_counter >= max_iterations:
                break
            iteration_counter += 1

            if iteration_counter == 4 or iteration_counter == 3:
                pass

            new_interpolation_point = self.new_interpolation_point(normalized_xy, t_values, x_curve, y_curve)

            #if 'suitable' values for a new interpolation point were found update control_points, control_weights, t_values, and bounds
            if new_interpolation_point is not None:
                control_points, control_weights, t_values, index_of_new_values = self.elevate_barycentric_curve(control_points, control_weights, t_values, new_interpolation_point)

                # bounds.insert(num_control_points + 2 * index_of_new_values, (control_points[2 * index_of_new_values + 1] / 1.05, control_points[2 * index_of_new_values + 1] * 1.05))
                # bounds.insert(num_control_points + 2 * index_of_new_values, (control_points[2*index_of_new_values]/1.05,control_points[2*index_of_new_values]*1.05))
                bounds.insert(num_control_points + 2 * index_of_new_values, y_bound)
                bounds.insert(num_control_points + 2 * index_of_new_values, x_bound)
                bounds.insert(0,weight_bounds)

                args_curve_error["t_values"] = t_values
                num_control_points += 1
            else:
                break

        x_curve, y_curve = self.scale_curve(x_curve, y_curve, 1 / sp.S(normalization_info["scale_factor"]))
        #curvature = self.calculate_curvature(x_curve, y_curve, t_values)[0] #calculates total curvature
        curvature = self.peak_curvature(x_curve, y_curve, t_values) #calculates peak curvature
        #curvature = curvature / normalization_info["scale_factor"]
        print("Normalized curvature: " + str(curvature[0]) + "\t" + "t_with_greatest_curvature: " + str(curvature[1]))
        return curvature[0]

    def new_interpolation_point(self, normalized_xy, t_values, x_curve, y_curve):
        small_val = 0.0001
        default_t = 0.5 - small_val
        is_left_to_right = normalized_xy[0][0] < normalized_xy[-1][0]
        t = sp.symbols("t", real=True)
        new_point = {"x":-1, "y":-1, "t":-1}

        # ensures x_curve & y_curve are copied appropriately
        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('j', 't'), local_dict={'t': t}), sp.parse_expr(str(y_curve).replace('j', 't'), local_dict={'t': t})

        #lambdify formulas for easy evaluation
        x_curve_lambda, y_curve_lambda = sp.lambdify(t, x_curve, modules="numpy"), sp.lambdify(t, y_curve, modules="numpy")
        dy_dx = sp.lambdify(t, y_curve.diff(t) / x_curve.diff(t), modules="numpy")

        #coord is the leftmost point of a segment of the polyline, slope is the slope of the line segment,
        # and length is how long the polyline segment is
        default_polyline_segment = {"startpoint":(-1,-1), "slope":0, "length":0, "endpoint":(-1,-1), "difference":0}
        polyline = []

        for i in range(len(normalized_xy) - 1):
            polyline.append(copy.deepcopy(default_polyline_segment))
            polyline[-1]["startpoint"] = normalized_xy[i]
            polyline[-1]["slope"] = (normalized_xy[i + 1][1] - normalized_xy[i][1]) / (normalized_xy[i + 1][0] - normalized_xy[i][0])
            polyline[-1]["length"] = sqeuclidean(normalized_xy[i], normalized_xy[i + 1])
            polyline[-1]["endpoint"] = normalized_xy[i + 1]


        def dot_product(ts:list[np.float64], point:list[tuple]):
            nonlocal x_curve_lambda, y_curve_lambda, dy_dx
            output = []

            for j in ts:
                if math.isnan(x_curve_lambda(j)) or math.isnan(y_curve_lambda(j)) or math.isnan(dy_dx(j)):
                    j += small_val

                curve_tangent = dy_dx(j)
                point_to_curve_vector = (point[0][0] - x_curve_lambda(j), point[0][1] - y_curve_lambda(j))
                output.append(abs(1 * point_to_curve_vector[0] + curve_tangent * point_to_curve_vector[1]))

            return output

        #finds the differences between the polyline and the current parametric curve and stores it in polylines difference attribute
        for i in range(len(polyline)):
            t0 = scipy_optimize.minimize(fun=dot_product, x0=default_t, args=[polyline[i]["startpoint"]], method="Nelder-Mead", bounds=[(small_val,1-small_val)]).x[0]
            t1 = scipy_optimize.minimize(fun=dot_product, x0=default_t, args=[polyline[i]["endpoint"]], method="Nelder-Mead", bounds=[(small_val,1-small_val)]).x[0]

            # list coords clockwise, required for shoelace formula
            if is_left_to_right:
                coords = [polyline[i]["startpoint"], polyline[i]["endpoint"], (x_curve_lambda(t1), y_curve_lambda(t1)),
                          (x_curve_lambda(t0), y_curve_lambda(t0))]
            else:
                coords = [polyline[i]["startpoint"], (x_curve_lambda(t0), y_curve_lambda(t0)),
                          (x_curve_lambda(t1), y_curve_lambda(t1)), polyline[i]["endpoint"]]

            #implementation of the shoelace formula
            area = 0
            for j in range(len(coords)):
                area += coords[j][1] * (coords[(j - 1) % 4][0] - coords[(j + 1) % 4][0])
            area = 0.5 * area

            polyline[i]["difference"] = area

        #where a region consists of consecutive segments from polyline that share the same difference's sign.
        default_region = {"total_difference":0, "included_segments":[0,0]}
        regions = []

        #populates regions
        total = 0
        for i in range(len(polyline)):
            #if polyline[i] is the last value in polyline
            if i == len(polyline) - 1:
                regions[-1]["total_difference"] = total
                regions[-1]["included_segments"][1] = i

            #if i > 0 and sign(polyline[i - 1]["difference"]) == sign(polyline[i]["difference"])
            elif i > 0 and polyline[i - 1]["difference"] / abs(polyline[i - 1]["difference"]) == polyline[i]["difference"] / abs(polyline[i]["difference"]):

                total += polyline[i]["difference"]

            #if sign(polyline[i - 1]["difference"]) != sign(polyline[i]["difference"])
            elif i > 0:
                regions[-1]["total_difference"] = total
                regions[-1]["included_segments"][1] = i - 1

                regions.append(copy.deepcopy(default_region))
                regions[-1]["included_segments"][0] = i

                total = polyline[i]["difference"]

            elif i == 0:
                regions.append(copy.deepcopy(default_region))
                total = polyline[i]["difference"]


        #find the largest absolute difference between polyline and the parametric curve
        abs_maximum, index = abs(regions[0]["total_difference"]), 0
        for i in range(1, len(regions)):
            if abs(regions[i]["total_difference"]) > abs_maximum:
                abs_maximum = abs(regions[i]["total_difference"])
                index = i

        if len(regions) >= 1 and regions[0]["total_difference"] != 0 and regions[0]["included_segments"] != [0,0]:
            index = int((regions[index]["included_segments"][0] + regions[index]["included_segments"][1]) / 2)
            closest_t = scipy_optimize.minimize(fun=dot_product, x0=default_t, args=[polyline[index]["startpoint"]],
                                                method="Nelder-Mead", bounds=[(small_val, 1 - small_val)]).x[0]

            new_point["t"], new_point["x"], new_point["y"] = closest_t, x_curve_lambda(closest_t), y_curve_lambda(closest_t)

        else:
            new_point = None

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
                interpolation_points.insert(2*i, new_interpolation_point_values["y"])
                interpolation_points.insert(2*i, new_interpolation_point_values["x"])
                
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

                if temp_weight < 0:
                    temp_weight = 0
                    sign = (-1) ** (len(weights) + new_interpolation_point_index)
                    for j in range(len(weights)):
                        temp_weight += sign * weights[j] / (t_values[j] - new_interpolation_point_values["t"])
                        sign = -sign

                new_weights.append(temp_weight)

        return interpolation_points, new_weights, new_t_values, new_interpolation_point_index

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

    def peak_curvature(self, x_curve:str, y_curve:str, t_values:list[np.float64]) -> tuple:
        """returns the peak of curvature over a parametric curve from 0 to 1, returns the calculated curvature and the error associated"""
        t = sp.symbols("t", real=True)
        small_val = 0.00001
        peak_curvature = 0

        for i in range(len(t_values)):
            t_values[i] = np.float64(t_values[i])

        #ensures x_curve & y_curve are copied appropriately
        x_curve, y_curve = sp.parse_expr(str(x_curve).replace('j','t'),local_dict={'t':t}), sp.parse_expr(str(y_curve).replace('j','t'),local_dict={'t':t})

        dx_dt, dy_dt = x_curve.diff(t), y_curve.diff(t) #finding the 1st derivatives
        ddx_dt, ddy_dt = sp.diff(dx_dt, t), sp.diff(dy_dt, t) #finding the 2nd derivatives

        numerator = abs(dx_dt * ddy_dt - dy_dt * ddx_dt)
        denominator = (dx_dt**2  + dy_dt**2)**sp.S(3/2)
        curvature = sp.lambdify(t, (numerator / denominator), modules="numpy")

        print("Curvature function: " + sp.latex(sp.S(numerator / denominator)))

        #finds the local maxima of the curvature function on the interval (0,1)
        list_of_t_values = np.linspace(small_val, 1, 997, endpoint=True)

        #prevents list_of_t_values from containing values in t_values, since values of t near any t_values may result in singularities
        i = 0
        while i < len(list_of_t_values):
            check = [math.isclose(list_of_t_values[i], j, abs_tol=small_val) for j in t_values]
            if True in check:
                list_of_t_values = np.delete(list_of_t_values, i)
            else:
                i += 1

        curvature_maxima = scipy_signal.argrelextrema(np.array([curvature(x) for x in list_of_t_values]), np.greater)[0]

        t_with_greatest_curvature = list_of_t_values[curvature_maxima[abs(np.array([list_of_t_values[i] for i in curvature_maxima]) - 0.5).argmin()]]
        peak_curvature = curvature(t_with_greatest_curvature)

        return peak_curvature, t_with_greatest_curvature

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
        interpolation_points = [1, 0, 12 / 13, 5 / 13, 0, 1]
        t_vals = [0, 1 / 3, 1]
        weights = [2, 13 / 6, 1 / 2]

        # temp = []
        # for i in range(int(len(interpolation_points)/2)):
        #     temp.append((interpolation_points[2*i],interpolation_points[2*i+1]))

        #x_curve, y_curve = test.barycentric_expression(temp, weights, t_vals)
        #print(sp.latex(sp.S(x_curve.replace('j', 't'))) + "\n" + sp.latex(sp.S(y_curve.replace('j', 't'))))
        new_point = {"x": 3 / 5, "y": 4 / 5, "t": 2 / 3}

        expected_interpolation_points = [1.0, 0.0, 12.0 / 13.0, 5.0 / 13.0, 3.0 / 5.0, 4.0 / 5.0, 0.0, 1.0]
        expected_t_vals = [0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0]
        expected_weights = [3.0, 13.0 / 2.0, 5.0, 1.4999999999999998]

        interpolation_points, weights, t_vals, trash = test.elevate_barycentric_curve(interpolation_points, weights, t_vals, new_point)
        interpolation_points, weights, t_vals = [float(x) for x in interpolation_points], [float(x) for x in weights], [float(x) for x in t_vals]

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
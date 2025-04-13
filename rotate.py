import math

########################################################################################
#    Author: Mark Dickinson - https://stackoverflow.com/users/270986/mark-dickinson
#    Date: 12/19/2015
#    Availability: https://stackoverflow.com/questions/34372480/rotate-point-about-another-point-in-degrees-python
########################################################################################

def rotate(origin: tuple, point: tuple, angle: float):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

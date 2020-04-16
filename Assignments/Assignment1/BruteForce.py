# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: ___________.py
#
# author: Peter Mackenzie-Helnwein
# created: 03/27/2020
#
# #######################################################

# import required libraries *****************************

import sys
from numpy import log, sqrt, linspace
import matplotlib.pyplot as plt


# defining functions and classes ************************
def F(u, v, EA=2100, d=0):

    # Calculate variables to be used for calculations
    l12 = (5.5 + u) ** 2 + (0.5 - v) ** 2
    l1 = sqrt(l12)
    L12 = 30.5
    L1 = sqrt(L12)

    l22 = (4 - u) ** 2 + (0.5 - v) ** 2
    l2 = sqrt(l22)
    L22 = 16.25
    L2 = sqrt(L22)

    # Calculate total force in x and y direction for given displacement
    if d == 0:
        fx1 = (5.5+u)/l1 * EA * 0.5 * log(l12/L12)
        fx2 = (u-4)/l2 * EA * 0.5 * log(l22/L22)
        fx = fx1+fx2

        fy1 = (0.5 - v) / l1 * EA * 0.5 * log(l12 / L12)
        fy2 = (0.5 - v) / l2 * EA * 0.5 * log(l22 / L22)
        fy = fy1 + fy2

        return fx, fy

    if d == 1:
        du = 1e-12
        dv = 1e-12

        fx, fy = F(u,v)

        fxdu1, fydu1 = F(u+du, v)

        fxdu = (fxdu1-fx) / du
        fydu = (fydu1-fy) / du


        fxdv1, fydv1 = F(u, v + dv)

        fxdv = (fxdv1 - fx) / dv
        fydv = (fydv1 - fy) / dv


        return fxdu, fydu, fxdv, fydv


def solver(x, y, P, tol=1e-12, m=1, data=False):



    sx = x
    sy = y

    xi = [x]
    yi = [y]

    Rx, Ry = F(sx, sy)

    countx = 0
    county = 0

    if m == 1:
        while abs(Ry - P) > tol and county < 100:
            countx = 0
            dFxdu, dFydu, dFxdv, dFydv = F(sx, sy, d=1)

            sy -= (Ry-P) / dFydv

            Rx, Ry = F(sx, sy)

            while abs(Rx) > tol and countx < 100:
                dFxdu, dFydu, dFxdv, dFydv = F(sx, sy, d=1)

                sx -= Rx / dFxdu

                Rx, Ry = F(sx, sy)
                countx += 1

            plt.plot(xi,yi)

            county += 1
    if m == 2:
        while (abs(Ry - P) > tol or abs(Rx) > tol) and county < 100:
            dFxdu, dFydu, dFxdv, dFydv = F(sx, sy, d=1)

            sy -= (Ry - P) / dFydv
            sx -= Rx / dFxdu

            xi.append(sx)
            yi.append(sy)

            Rx, Ry = F(sx, sy)


            print('{:3d}: x={:16.12f}  Fy(u,v)={:16.12e}  Fx(u,v)={:16.12e}'.format(county, sy, Ry-P, Rx))
            county += 1


    if data:
        plt.plot(xi, yi)
        plt.show()

    return sx, sy












# defining main execution procedure
def main():
    pass


# main execution ****************************************

if __name__ == "__main__":
    print(F(0, 0.1))
    print(F(0, 0.1, d=1))
    '''
    p = linspace(0, -0.9, 100)
    xi = []
    yi = []

    for i in p:
        x, y = solver(0, 0, i, m=2)

        xi.append(x)
        yi.append(y)

    plt.plot(xi, yi)
    plt.show()
    '''

    u, v = solver(0, 0, -0.5, tol=1e-12, m=2, data=True)

    print([u, v])


    sys.exit(0)
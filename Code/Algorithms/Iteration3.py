# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: Iteration1.py
#
# author: Peter Mackenzie-Helnwein
# created: 03/27/2020
#
# #######################################################

'''
This file is to demonstrate the iterative solution of
a nonlinear equation

Version 3: Newton method
'''

# import required libraries *****************************

import sys
import numpy as np

# import everything from my plotter class
from Plotter import *

# defining functions and classes ************************

def F(x, d=0):
    if d == 1:
        ans = 2*x
    else:
        ans = x*x - 10.
    return ans


def solver3(startX=0.0, TOL=1.0e-12):
    '''
    This solver uses Newton's method to find the root of F(x)
    starting from a given estimate.

    The solver uses a default tolerance of 10^(-12).
    You may overwrite the tolerance by specifying, e.g.,
    TOL=1.0e-5 when calling this solver.
    '''

    # record starting point
    interval = [startX]

    # iterating for the root

    s = startX
    ys  = F(s)  # anything larger than TOL will do
    xi = [s]
    yi = [ys]

    cnt = 0         # set counter for emergency exit in case no solution is found

    print('{:3d}: x={:16.12f}  F(x)={:16.12e}'.format(cnt,s,ys))

    while abs(ys) > TOL and cnt < 100:
        dydx = F(s, 1)

        s -= ys / dydx
        ys = F(s)

        xi.append(s)
        yi.append(ys)

        cnt += 1

        print('{:3d}: x={:16.12f}  F(x)={:16.12e}'.format(cnt,s,ys))

    # plotting the nonlinear function over that interval
    x = np.linspace(np.min(xi), np.max(xi), 1000)
    y = F(x)

    # collect iteration data
    answer = {
        'x':x,
        'y':y,
        'xi':xi,
        'yi':yi,
        'label':'Newton Algorithm',
        'save':'Iteration3.png'
    }

    # return the answer
    return answer


def main():
    ans = solver3(1.)
    Plotter(**ans, tangent=True, path=False)


# main execution ****************************************

main()
sys.exit(0)

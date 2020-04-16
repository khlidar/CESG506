'''
Created by: Kristinn Hlíðar Grétarsson
Date:       04.12.2019

Description:

Functions:

'''

# Imports
import sys
from SolverV2 import *
# Functions and classes

# defining main execution procedure
def main():
    X1 = array([0., 0.])
    X2 = array([5.5, 0.5])
    X3 = array([9.5, 0.])

    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}

    ele1 = TrussElement((X1, X2), para)
    ele2 = TrussElement((X2, X3), para)

    p = array([0., -0.981713439866848])
    l = [0, 0.25, 0.5, 0.75, 0.99, 0.999]


    solve(ele1, ele2, p, lf=l, data=True)
    pass


# main execution ****************************************
if __name__ == "__main__":
    main()
    sys.exit(0)
# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: ___________.py
#
# author: Kristinn Hlíðar
# created: 2020-04-13
#
# #######################################################

# import required libraries *****************************

import sys
from truss import *
from numpy import array, linalg, asarray, transpose

# defining functions and classes ************************
def solve(ele1:TrussElement, ele2:TrussElement, P, lf=[0, 1], tol=1e-12):
    ul = [] # keeps track of disp. at different load steps

    for lfp in lf:
        K1 = ele1.getKe()
        K2 = ele2.getKe()
        Ktot = K1[1, 1] + K2[0, 0]

        R1 = ele1.getFeAsMatrix()
        R2 = ele2.getFeAsMatrix()

        Rtot = R1[2:4] + R2[0:2]
        Rtot = asarray(Rtot)


        Ptot = lfp * P
        u = linalg.solve(Ktot, Ptot - Rtot)
        u = transpose(u)
        u = u.flatten()

        nulV = array([0., 0.])

        ele1.setDisp((nulV, array([u[0], u[1]])))
        ele2.setDisp((array([u[0], u[1]]), nulV))

        R1 = ele1.getFeAsMatrix()
        R2 = ele2.getFeAsMatrix()

        Rtot = R1[2:4] + R2[0:2]
        Rtot = asarray(Rtot)

        count = 0
        while (abs(Rtot[0] - Ptot[0]) > tol or abs(Rtot[1] - Ptot[1]) > tol) and count < 100:
            # Get stiffness of element 1 and two
            K1 = ele1.getKe()
            K2 = ele2.getKe()
            # Assemble stiffness at free node (assumes node 2 for ele1 and node 1 for ele2)
            Ktot = K1[1, 1] + K2[0, 0]

            # Displacement step
            du = linalg.solve(Ktot, Ptot - Rtot)
            du = transpose(du)
            du = du.flatten()

            # Add displacement step to total displacement
            u += du

            # Add displacement to elements
            ele1.setDisp((nulV, array([u[0], u[1]])))
            ele2.setDisp((array([u[0], u[1]]), nulV))

            # Get forces in elements
            R1 = ele1.getFeAsMatrix()
            R2 = ele2.getFeAsMatrix()
            # Get total force in node by adding together node 2 element 1 and node 1 element 2
            Rtot = R1[2:4] + R2[0:2]
            Rtot = asarray(Rtot)

            # Add to counter
            count += 1

        ul.append(u)

    return u, ul






# defining main execution procedure

def main():
    X1 = array([0., 0.])
    X2 = array([5.5, 0.5])
    X3 = array([9.5, 0.])

    para = {'E': 2100.0,
            'A': 1.0,
            'nu':0.0}
    
    ele1 = TrussElement((X1, X2), para)
    ele2 = TrussElement((X2, X3), para)

    p = array([[0.], [-0.981]])
    l = [0, 0.25, 0.5, 0.75, 0.99, 0.999]

    print(solve(ele1, ele2, p, lf=l))
    pass


# main execution ****************************************

if __name__ == "__main__":
    main()
    sys.exit(0)
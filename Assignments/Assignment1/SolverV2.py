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
from numpy import array, linalg, asarray, transpose, linspace

# defining functions and classes ************************
def solve(ele1:TrussElement, ele2:TrussElement, P, lf=[0, 1], tol=1e-12, data=False):
    extraData = []
    ul = [] # keeps track of disp. at different load steps

    u = array([0., 0.])

    # Create a zero array, used when assigning displacements
    nulV = array([0., 0.])

    for lfp in lf:
        # Get force to solve for
        Ptot = lfp * P

        # Get updated reaction forces
        Rtot = force(ele1, ele2, u)

        count = 0
        if data:
            print(f'Load step = {lfp}')
            error = Ptot - Rtot
            error = linalg.norm(error)
            print(f'{count}, {error:.6e}')
            #print(f'#{count} u={u[0]:.4e} v={u[1]:.4e} Rx={(Rtot - Ptot)[0]:10.6e} Ry={(Rtot - Ptot)[1]:10.6e} ')

        while (abs(Rtot[0] - Ptot[0]) > tol or abs(Rtot[1] - Ptot[1]) > tol) and count < 100:
            # Get stiffness matrix for free degree of freedom
            Ktot = stiffness(ele1, ele2, u)

            # Displacement step
            du = linalg.solve(Ktot, Ptot - Rtot)

            # Add displacement step to total displacement
            u += du

            # Get forces in elements
            Rtot = force(ele1, ele2, u)

            # Add to counter
            count += 1
            if data:
                # Calculate error
                error = Ptot - Rtot
                error = linalg.norm(error)
                print(f'{count}, {error:.6e}')
                #print(f'#{count} u={u[0]:.4e} v={u[1]:.4e} Rx={(Rtot - Ptot)[0]:10.6e} Ry={(Rtot - Ptot)[1]:10.6e} ')
                extraData.append([count, lfp, error])

        if data:
            print("")
        ul.append(u)

    return u, extraData

def force(ele1:TrussElement, ele2:TrussElement, u=array([0., 0.])):
    # Create a zero array, used when assigning displacements
    nulV = array([0., 0.])

    ele1.setDisp((nulV, u))
    ele2.setDisp((u, nulV))

    # Get updated reaction forces
    R1 = ele1.getFe()
    R2 = ele2.getFe()
    Rtot = R1[1] + R2[0]

    return Rtot

def stiffness(ele1:TrussElement, ele2:TrussElement, u=array([0., 0.])):
    # Create a zero array, used when assigning displacements
    nulV = array([0., 0.])

    ele1.setDisp((nulV, u))
    ele2.setDisp((u, nulV))

    # Get stiffness of element 1 and two
    K1 = ele1.getKe()
    K2 = ele2.getKe()
    # Assemble stiffness at free node (assumes node 2 for ele1 and node 1 for ele2)
    Ktot = K1[1, 1] + K2[0, 0]

    return Ktot



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
    p = array([0., -0.981713439866848])
    l = [0, 0.25, 0.5, 0.75, 0.99, 0.999]

    #print(solve(ele1, ele2, p, lf=l, data=True))

    test = solve(ele1, ele2, p, lf=l, data=True)

    '''
    # Finding P.cr
    pcr = 0

    X1 = array([0., 0.])
    X2 = array([5.5, 0.5])
    X3 = array([9.5, 0.])

    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}

    ele1 = TrussElement((X1, X2), para)
    ele2 = TrussElement((X2, X3), para)

    u = array([0., -0.2])

    u0 = linspace(-0.21, -0.22, 200)
    tol = 1e-4
    for i in u0:
        u = array([0., i])

        Rtot = force(ele1, ele2, u)
        rx = Rtot[0]

        count = 0
        while abs(rx) > tol:
            Ksys = stiffness(ele1, ele2, u)

            K11 = Ksys[0, 0]
            K12 = Ksys[0, 1]

            ux = rx / K11

            u[0] -= ux

            Rtot = force(ele1, ele2, u)
            rx = Rtot[0]

            count += 1

            if count > 20:
                print("not working")
                break

        if Rtot[1] < pcr:
            pcr = Rtot[1]

        print(f'{u},{Rtot[1]}')

    print(pcr)
    '''

    pass


# main execution ****************************************

if __name__ == "__main__":
    main()
    sys.exit(0)
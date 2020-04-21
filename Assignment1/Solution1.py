# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: Solution1.py
#
#    solutions for assignment #1
#
# author: Peter Mackenzie-Helnwein
# contributor: William(Bill) and Tatsu (P_critical for Problem 1-2)
# created: 04/10/2020
# updated: 04/14/2020
#
# #######################################################

# import required libraries *****************************

import sys
import numpy as np

from truss import *

import matplotlib.pyplot as plt

# defining functions and classes ************************

def getForce_1_1(vertical_disp):
    """
        get a vector of reaction forces for Problem 1-1

        the answer is organized as follows:

        idx   strain     equilibrium
        ----------------------------
        0     linear     undeformed
        1     Green      undeformed
        2     Almansi    undeformed
        3     Henkey     undeformed
        ----------------------------
        4     linear     deformed
        5     Green      deformed
        6     Almansi    deformed
        7     Henkey     deformed
        ----------------------------

        I placed the origin of my coordinate system underneath the sliding (right) support
        at the height of the lower (left) support.

        x points to the right, y points up.
        u is positive if the right support moves down.
        P is positive pointing downward.
        """

    # Kinematics
    # ==========

    # Length, axial vectors

    XI = np.array([-5.5, 0])
    XJ = np.array([0, 0.5])
    ui = np.array([0, 0])
    uj = np.array([0, -vertical_disp])

    Lvec = XJ - XI

    xi = XI + ui
    xj = XJ + uj

    lvec = xj - xi

    Len2 = Lvec @ Lvec

    L = np.sqrt(Len2)

    len2 = lvec @ lvec

    len = np.sqrt(len2)

    Nvec = Lvec / L

    nvec = lvec / len

    # Strain & Stretch

    lambda2 = len2 / Len2

    # Traditional definition (similar to linear theory, but nonlinear)

    Epsilon1 = np.sqrt(lambda2) - 1

    # Green-Lagrange strain

    Epsilon2 = 0.5 * (lambda2 - 1.)

    # Almansi strain

    Epsilon3 = 0.5 * (1. - 1./lambda2)

    # Henkey or Natural strain

    Epsilon4 = 0.5 * np.log(lambda2)

    strain = np.array([Epsilon1, Epsilon2, Epsilon3, Epsilon4])

    # Constitutive
    # ============

    # We assume that the cross section area does not change with deformation.
    # This is not always a good assumption but it is viable with our very small strains.

    EA = 2100

    F = EA * strain

    # Equilibrium
    # ===========

    # On the undeformed system

    P0 = -F * Nvec[1]  # the y-component of N

    # On the deformed system

    Pn = -F  * nvec[1]  # the y-component of n

    return (strain, np.hstack((P0, Pn)))

def NewtonRaphson():
    pass

def Problem_1_1():
    """
    Problem 1-1

    I placed the origin of my coordinate system underneath the sliding (right) support
    at the height of the lower (left) support.

    x points to the right, y points up.
    u is positive if the right support moves down.
    P is positive pointing downward.
    """

    # compute results for -0.1 h <= u <= 2.1 h
    # ========================================

    h = 0.5  # height of the truss

    forces = []
    strains = []

    u =np.linspace(-0.1*h, 2.1*h, 50)
    #u =np.linspace(-0.*h, 2.5*h, 50)

    for disp in u:
        (eps, P) = getForce_1_1(disp)
        forces.append(P)
        strains.append(eps)

    forces  = np.array(forces)
    strains = np.array(strains)

    # plotting results
    # ================

    fig, (ax1, ax2, ax3) = plt.subplots(1,3)

    # set colors for plots
    type = ['r-','g-','b-','c-']
    labels = ["Traditional", "Green-Lagrange", "Almansi", "Henkey"]

    # plotting strains
    # ------------------------------------------------
    for i in range(4):
        ax1.plot(u, strains[:,i],type[i],label=labels[i])
    ax1.set_title('kinematics')
    ax1.set_xlabel('displacement, $u$')
    ax1.set_ylabel('strain, $\epsilon$')
    ax1.legend()
    ax1.grid(True)

    # Reality check: the horizontal configuration should have a strain in the range of (5.5 - L)/L

    # This comparison shows that for small strains, all options are viable and
    # yield nearly identical results (GREAT!)

    # Plots using equilibrium on the undeformed system
    # ------------------------------------------------
    for i in range(4):
        ax2.plot(u, forces[:,i],type[i],label=labels[i])
    ax2.set_title('undeformed system')
    ax2.set_xlabel('displacement, $u$')
    ax2.set_ylabel('force, $P$')
    ax2.legend()
    ax2.grid(True)

    # Note how all options predict zero force at u=2H, which completely relaxes the bar.
    # However, there shouldn't be a vertical force as u=H levels the bar.

    # Plots using equilibrium on the deformed system
    # ------------------------------------------------
    for i in range(4):
        ax3.plot(u, forces[:,i+4],type[i],label=labels[i])
    ax3.set_title('deformed system')
    ax3.set_xlabel('displacement, $u$')
    ax3.set_ylabel('force, $P$')
    ax3.legend()
    ax3.grid(True)

    # Note how all options predict zero force at u=H and at u=2H.

    # This problem is one where displacements are very large but the associated strain
    # remains quite small.  As a consequence, all of the discusses strain measures
    # yield suitable and extremely similar results.
    # The confidence level in the solution to this problem is thus quite high.

    # export figure to file
    fig.set_size_inches(12, 4, forward=True)
    plt.savefig('problem1-1_solution.png', dpi=100)


class Problem_1_2(object):
    """
        Problem 1-2

        I placed the origin of my coordinate system underneath the sliding (right) support
        at the height of the lower (left) support.

        x ... points to the right, y points up.
        u ... is positive it moves in the positive x-direction.
        v ... is positive it moves in the positive y-direction.
        Px .. is positive pointing in the positive x-direction.
        Py .. is positive pointing in the positive y-direction.
    """

    TOL = 5.0e-12

    def __init__(self):

        # setup
        # ==========

        # nodes

        X1 = np.array([-5.5, 0])
        X2 = np.array([0, 0.5])
        X3 = np.array([4.0, 0.0])

        # elements (not needed for the brute force implementation)

        # ... I define both elements such that node i is the free node
        #     This simplifies assembly of [Kss] and {Fs}
        # ... since only the product EA matters, set E=2100 and A=1
        self.elem1 = TrussElement((X2, X1), params={'E': 2100.0, 'A': 1.0})
        self.elem2 = TrussElement((X2, X3), params={'E': 2100.0, 'A': 1.0})

        # set initial solution for the free node

        # sol is a list of equilibrium states
        # each equilibrium state is representated as dictionary with keys
        #    lambda ... the dimensionless load factor
        #    U ........ the displacement vector of the free node (as np.array)

        self.sol = [{'lambda': 0.0, 'U': np.array([0.0, 0.0])}]

        # define reference load vector on free node
        #self.Pref = np.array([0.0, -0.990])      # my quick and dirty critical load
        self.Pref = np.array([0.0, -0.98171345])  # thanks to William and Tatsu for the more accurate value

    def __str__(self):
        return "this class solves problem 2 from Assignment #1"

    def newtonRaphson(self, loadfactor, U0=np.array([0.0,0.0])):
        # initialize
        U = U0.copy()

        # check for equilibrium (this returns norm of R!)
        normR = self.getResiduum(loadfactor, U)
        normR0 = np.max((normR, 0.1))

        # provide useful information at initialization
        print("# {:3d}: U=({:12.6g},{:12.6g}) err={:16.6e}".format(0, U[0], U[1], normR/normR0))

        # perform iteration loop
        cnt = 1   # counter for emergency exit
        while normR/normR0 > self.TOL:
            # update U (single iteration step)
            U += np.linalg.solve(self.KT, self.R)

            # check for equilibrium (this returns norm of R!)
            normR = self.getResiduum(loadfactor, U)             # the object oriented approach
            #normR = self.getResiduumBruteForce(loadfactor, U)   # the brute-force approach

            # provide useful information from iteration steps
            print("# {:3d}: U=({:12.6g},{:12.6g}) err={:16.6e}".format(cnt, U[0], U[1], normR/normR0))

            # check if we exceed the max number of steps == emergency exit
            cnt += 1
            if cnt>10:
                print("*** warning: iteration failed to converge")
                break

        # add empty line
        print()

        # return result
        return U

    def getResiduumBruteForce(self, lf, disp):

        # external load

        self.R  = lf * self.Pref
        self.KT = np.zeros((2,2))

        # element a

        Lvec = np.array([5.5, 0.5])
        Ui = np.zeros_like(disp)
        Uj = disp

        # ... independent of the element

        Le2 = np.dot(Lvec, Lvec)

        lene = Lvec + Uj - Ui
        le2 = np.dot(lene, lene)
        ne = lene / np.sqrt(le2)

        eps = 0.5 * np.log(le2 / Le2)

        EA = 2100
        fe = EA * eps

        self.R -= fe * ne

        self.KT += EA / np.sqrt(le2) * np.outer(ne, ne)
        self.KT += fe / np.sqrt(le2) * (np.identity(2) - np.outer(ne, ne))

        # element b

        Lvec = np.array([-4.0, 0.5])
        Ui = np.zeros_like(disp)
        Uj = disp

        # ... independent of the element

        Le2 = np.dot(Lvec, Lvec)

        lene = Lvec + Uj - Ui
        le2 = np.dot(lene, lene)
        ne = lene / np.sqrt(le2)

        eps = 0.5 * np.log(le2 / Le2)

        EA = 2100
        fe = EA * eps

        self.R -= fe * ne

        self.KT += EA / np.sqrt(le2) * np.outer(ne, ne)
        self.KT += fe / np.sqrt(le2) * (np.identity(2) - np.outer(ne, ne))

        return np.sqrt(self.R @ self.R)

    def getResiduum(self, lf, disp):

        self.elem1.setDisp((disp, np.zeros_like(disp)))
        self.elem2.setDisp((disp, np.zeros_like(disp)))

        self.R = lf * self.Pref
        self.R -= self.elem1.getFe()[0]
        self.R -= self.elem2.getFe()[0]

        self.KT = np.zeros((2, 2))
        self.KT += self.elem1.getKe()[0, 0]
        self.KT += self.elem2.getKe()[0, 0]

        return np.sqrt(self.R @ self.R)

    def laodcontrol(self, loadfactors):# start value for Newton-Raphson iteration

        for lf in loadfactors:

            print("*** new load factor: {} ***\n".format(lf))

            # set start value
            U0 = self.sol[-1]['U']

            # perform iterative solution algorithm
            U = self.newtonRaphson(lf, U0)

            # add found solution to the solution list (equilibrium path)
            self.sol.append({'lambda': lf, 'U': U})

    def report(self):
        # this would be a way to generate a nice report

        # data collection arrays for plotting
        lf = []
        Ux = []
        Uy = []

        # print equilibrium path information
        print("  --- lf ---   ----- Ux -----   ----- Uy -----")
        for state in self.sol:

            lf.append(state['lambda'])
            Ux.append(state['U'][0])
            Uy.append(state['U'][1])

            print("{:12.5f},{:16.10f},{:16.10f}".format(state['lambda'],*list(state['U'])))

        # == plotting results ==

        fig, (ax1, ax2) = plt.subplots(1, 2)

        # plotting the equilibrium path
        ax1.plot(Ux,lf,'mo-',label="$U_x\equiv u$ (pos $\\rightarrow$)")
        ax1.plot(Uy,lf,'co-',label="$U_y\equiv v$ (pos $\\uparrow$)")
        ax1.grid(True)
        ax1.set_ylabel("load factor, $\lambda$")
        ax1.set_xlabel("displacement")
        ax1.legend()
        ax1.set_title("Load-Displacement curves")

        # plotting the equilibrium path
        ax2.plot(Ux,Uy,'ro-',label="trace of loaded node")
        ax2.grid(True)
        ax2.set_xlabel("displacement component $U_x\equiv u$")
        ax2.set_ylabel("displacement component $U_y\equiv v$")

        for x in zip(lf, Ux, Uy):
            ax2.annotate("{:8.3f}".format(x[0]), (x[1], x[2]))

        #ax2.legend()
        ax2.set_title("Tracing the loaded node")

        # export figure to file
        fig.set_size_inches(12, 6, forward=True)
        plt.savefig('problem1-2_solution.png',dpi=100)



# defining main execution procedure

def main():

    # solve problem 1-2
    # ======================
    Problem_1_1()

    # solve problem 1-2
    # ======================

    prob2 = Problem_1_2()
    print(prob2)
    # set load factors of interest
    lfs = [0,0.25,0.5,0.75,0.9,0.97,0.99,0.999]

    # compute for those load factors
    prob2.laodcontrol(lfs)

    # create a report
    prob2.report()


# main execution ****************************************

if __name__ == "__main__":
    main()
    sys.exit(0)

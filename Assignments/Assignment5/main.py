'''
Created by: Kristinn Hlíðar Grétarsson
Date:       2020.04.28

Description:
    Solution for assignment 3 CESG506

Functions:

'''

# Imports
from Structure import *
from numpy import linalg, log
import matplotlib.pyplot as plt

# Functions and classes
def problem3_1():
    # Create a structure for the problem
    p3 = struct()

    # Add nodes to structure
    p3.addNode(1, array([0, 0]), array([0., 0.]))
    p3.addNode(2, array([0, 1]), array([5.5, 0.5]))

    # Add element to structure
    para = {'E': 2100.0,
                'A': 1.0,
                'nu': 0.0}

    p3.addEle(1, 1, 2, TrussElement(params=para))
    p3.updateElements()

    # Load the structure
    P = 0.3
    p3.addLoad(1, 2, array([0., -P]))

    # Solve the structure
    alfa = 1.0
    upath, load, s_hist, R_hist = p3.solveArc(load=2, la=0.25, data=True, alfa=alfa)


    # Solution from assignment 1
    EA = 2100
    h = 0.5
    u = linspace(0, -1.2, 100)
    L = sqrt(30.5)
    l = sqrt(L**2 + u**2 + u)

    P1 = -(EA * log(l/L) * (h + u)/l)



    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    ax1.scatter(upath[:, 0], P * load, color='b', label='Arc length code')
    ax1.plot(u, P1, label='Assignment 1 solution')
    plt.title(f'Arc lenght check, alfa={alfa}')
    ax1.legend()

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)

    ax1.semilogy(range(len(R_hist)), R_hist, label='|R|')
    plt.title(f'|R| - iteration steps, alfa={alfa}')
    ax1.set_xlabel('Iteration steps')
    ax1.set_ylabel('|R|')
    ax1.legend()

    fig3 = plt.figure()
    ax1 = fig3.add_subplot(111)

    ax1.semilogy(s_hist, R_hist, label='|R|')
    ax1.set_xlabel('arc length')
    ax1.set_ylabel('|R|')
    plt.title(f'|R| - arc length, alfa={alfa}')
    ax1.legend()


    plt.show()


def problem3_2():
    p3 = struct()

    p3.addNode(1, array([0, 0]), array([0., 0.]))
    p3.addNode(2, array([1, 1]), array([5.5, 0.5]))
    p3.addNode(3, array([0, 0]), array([9.5, 0.]))

    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}

    p3.addEle(1, 1, 2, TrussElement(params=para))
    p3.addEle(1, 2, 3, TrussElement(params=para))

    p3.updateElements()

    p3.addLoad(1, 2, array([0., -1]))

    upath, load = p3.solveArc(load=2, la=0.2)



    # Assignment 2 solution
    p1 = struct()
    p1.addNode(1, array([0, 0]), array([0.0, 0.0]))
    p1.addNode(2, array([1, 1]), array([5.5, 0.5]))
    p1.addNode(3, array([0, 0]), array([9.5, 0.0]))

    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}

    p1.addEle(1, 1, 2, TrussElement(params=para))
    p1.addEle(2, 2, 3, TrussElement(params=para))

    p1.updateElements()

    p = p1.solveDis(array([2, array([0, 1])]), g=0.01)

    load2, upath2 = p1.displacementpath(array([2, array([0, -1])]), [0, -1.2], steps=100)
    load2 = array(load2)





    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(upath2[:, 1], load2, color='r', label='Displacement')
    ax1.scatter(upath[:, 1], load, color='b', label='Arc length')
    plt.title('Load vs vertical displacement')
    ax1.set_xlabel('Vertical disp. [m]')
    ax1.set_ylabel('Load factor')
    ax1.legend()

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    ax1.scatter(upath2[:, 0], load2, color='r', label='Displacement')
    ax1.scatter(upath[:, 0], load, color='b', label='Arc length')
    plt.title('Load vs horizontal displacement')
    ax1.set_xlabel('Horizontal disp. [m]')
    ax1.set_ylabel('Load factor')
    ax1.legend()

    plt.show()



def problem3_3():
    height = 5
    stories = 11

    width = height/20

    h = height/stories

    p3 = struct()

    para2000 = {'E': 2000.0,
                'A': 1.0,
               'nu': 0.0}

    para5000 = {'E': 5000.0,
                'A': 1.0,
                'nu': 0.0}

    j = 1
    for i in range(1, stories + 2):
        if i == 1:
            p3.addNode((i-1)*2+1, array([0, 0]), array([0., 0.]))
            p3.addNode((i-1)*2+2, array([0, 0]), array([width, 0.]))
        else:
            p3.addNode((i-1)*2+1, array([1, 1]), array([0., h*(i-1)]))
            p3.addNode((i-1)*2+2, array([1, 1]), array([width, h*(i-1)]))

            # Add vertical element
            p3.addEle(j, (i - 2) * 2 + 1, (i - 1) * 2 + 1, TrussElement(params=para2000))
            j += 1
            p3.addEle(j, (i - 2) * 2 + 2, (i - 1) * 2 + 2, TrussElement(params=para2000))
            j += 1
            # Add diagonal element
            p3.addEle(j, (i - 2) * 2 + 1, (i - 1) * 2 + 2, TrussElement(params=para5000))
            j += 1
            # Add horizontal element
            p3.addEle(j, (i - 1) * 2 + 1, (i - 1) * 2 + 2, TrussElement(params=para5000))
            j += 1

    P = 6.1685
    p3.addLoad(1, 23, array([0., -P/2]))
    p3.addLoad(2, 24, array([0., -P/2]))

    p3.getNodes(True)
    p3.getElements(True)

    p3.updateElements()
    p3.calcK()

    p3.plotStructure(mark=True)

    upath, load = p3.solveArc(s_lim=20, la=.2, load=10.0, nodeDis=[24, array([0., -5])], alfa=0.2)

    p3.plotStructure(True)


    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(upath[:, 40], load, color='b', label='Node$23_x$')
    ax1.scatter(upath[:, 42], load, color='r', label='Node$24_x$')
    plt.title(f'LF-disp. P={P} kN')
    ax1.set_xlabel('Horizontal displacement [m]')
    ax1.set_ylabel('Load factor')
    ax1.legend()

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(upath[:, 41], load, color='b', label='Node$23_y$')
    ax1.scatter(upath[:, 43], load, color='r', label='Node$24_y$')
    plt.title(f'LF-disp. P={P} kN')
    ax1.set_xlabel('Vertical displacement [m]')
    ax1.set_ylabel('Load factor')
    ax1.legend()

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    ax1.scatter(upath[:, 40], upath[:, 41]+5, color='b', label='Node 23')
    ax1.scatter(upath[:, 42]+.25, upath[:, 43]+5, color='r', label='Node 24')
    plt.title(f'Trace plot')
    ax1.set_xlabel('Horizontal location [m]')
    ax1.set_ylabel('Vertical location [m]')
    ax1.set(xlim=(0, 5.5), ylim=(0, 5.5))
    ax1.legend()


    plt.show()



# main execution ****************************************
if __name__ == "__main__":
    problem3_1()
    problem3_2()
    problem3_3()
    sys.exit(0)
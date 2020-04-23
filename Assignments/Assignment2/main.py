'''
Created by: Kristinn Hlíðar Grétarsson
Date:       2020.04.22

Description:
    FEM structure program that uses FEM elements from Peter to construct and calculate non-linear problems

Functions:

'''

# Imports
from Structure import *
from numpy import linalg
import matplotlib.pyplot as plt

# Functions and classes
def problem2_1():
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

    load, upath = p1.displacementpath(array([2, array([0, -1])]), [0, -1.2], steps=100)
    load = array(load)

    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    ax2.scatter(upath[:, 0], load, color='b', label='$u_x$')
    ax2.scatter(upath[:, 1], load, color='r', label='$u_y$')
    ax2.legend()
    ax2.set_xlabel('Displacement [m]')
    ax2.set_ylabel('Load factor')
    plt.title('Load factor vs displacement')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.scatter(upath[:, 0], upath[:, 1])
    ax2.set_xlabel('Horizontal displacement [m]')
    ax2.set_ylabel('Vertical displacement [m]')
    plt.title('Equilibrium path')


    X = linspace(-0.01, 0.01, 100)
    Y = linspace(0.2, -1.2, 100)
    Z = zeros((len(X), len(Y)))

    for i in range(len(X)):
        for j in range(len(Y)):
            p1.displaceEle([2, 3], array([X[i], Y[j]]))
            Z[j][i] = linalg.norm(p1.rtot[[2, 3]])

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    ax.contour(X, Y, Z)
    ax.scatter(upath[:, 0], upath[:, 1])
    plt.title('Total force for given displacement, contour')
    plt.show()


def problem2_2():
    # Create structure
    p2 = struct()

    # Add nodes to structure
    p2.addNode(1, array([0, 0, 0]), array([0.0, 5., 0.]))
    p2.addNode(2, array([0, 0, 0]), array([9.5, 5., 0.]))
    p2.addNode(3, array([0, 0, 0]), array([0.0, 0., 0.]))
    p2.addNode(4, array([0, 0, 0]), array([9.5, 0., 0.]))
    p2.addNode(5, array([1, 1, 1]), array([5.5, 3.75, 0.5]))
    p2.addNode(6, array([1, 1, 1]), array([5.5, 2.50, 0.5]))

    # Create dictionary for element parameters
    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}

    # Add elements to structure
    p2.addEle(1, 1, 5, TrussElement(params=para))
    p2.addEle(2, 1, 6, TrussElement(params=para))
    p2.addEle(3, 2, 5, TrussElement(params=para))
    p2.addEle(4, 3, 6, TrussElement(params=para))
    p2.addEle(5, 4, 6, TrussElement(params=para))
    p2.addEle(6, 4, 5, TrussElement(params=para))
    p2.addEle(7, 5, 6, TrussElement(params=para))

    # Update all elements
    p2.updateElements()

    # Update stiffness and internal forces
    p2.calcK()
    p2.partitionK()
    p2.calcR()

    # Get load and load path for enforced displacement of node 5
    load, upath = p2.displacementpath(array([5, array([0, 0, 1])]), [0, -1.2], steps=200)
    load = array(load)

    # Plot 3D displacement of node 5 and 6 through space
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.scatter(upath[:, 0]+5.5, upath[:, 1]+3.75, upath[:, 2]+0.5, color='b', label='$node 5$')
    ax1.scatter(upath[:, 3]+5.5, upath[:, 4]+2.5, upath[:, 5]+0.5, color='r', label='$node 6$')
    ax1.legend()
    plt.title('Location of nodes 5 & 6 through snap through')

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(121)
    ax2.scatter(upath[:, 0], load, color='b', marker='+', label='$u5_x$')
    ax2.scatter(upath[:, 1], load, color='r', marker='+', label='$u5_y$')
    ax2.scatter(upath[:, 3], load, color='b', label='$u6_x$')
    ax2.scatter(upath[:, 4], load, color='r', label='$u6_y$')
    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Displacement [m]')
    ax2.legend()

    ax2 = fig2.add_subplot(122)
    ax2.scatter(upath[:, 2], load, color='g', marker='+', label='$u5_z$')
    ax2.scatter(upath[:, 5], load, color='g', label='$u6_z$')
    ax2.set_ylabel('Load factor')
    ax2.set_xlabel('Displacement [m]')
    ax2.legend()
    plt.show()







# main execution ****************************************
if __name__ == "__main__":
    #problem2_1()
    problem2_2()
    sys.exit(0)
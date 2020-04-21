# #######################################################
#
# CESG 506 - Analysis of Nonlinear Structures
#
# file: Element.py
#
# author: Peter Mackenzie-Helnwein
# created: 04/05/2020
#
# #######################################################

# import required libraries *****************************

import sys
import copy
import numpy as np

# defining functions and classes ************************

class Element(object):
    """
    defining a master class for finite elements

    variables:
        self.X = tuple of nodal position vectors (as np.array)
        self.params = {'E':1.0, 'nu':0.0}
        self.nnode ......... number of nodes per element
        self.ndof .......... number of degrees of freedom per node
        self.U ............. array of nodal displacement vectors
        self.force ......... internal forces at nodes
        self.stiffness ..... the stiffness matrix

    methods:
        __init__(self, X, params)
        setDisp(self, U) ..... U is an array of nodal displacement vectors
        getFe(self) .......... return the internal force vector as array of nodal vectors
        getKe(self) .......... return the stiffness matrix as array of nodal matrices
        getFeAsMatrix(self) .. return the internal force vector as nx1 matrix
        getKeAsMatrix(self) .. return the stiffness matrix as nxn matrix
        def init(self) ....... element specific initialization steps
        compute(self) ........ does the actual computation
    """

    def __init__(self, X=(np.zeros(2),np.zeros(2)), params={'E':1.0, 'nu':0.0}):

        # check dimension
        self.nnode = len(X)
        ndof = 0
        for point in X:
            if ndof:
                pass
            else:
                ndof = len(point)
        self.ndof = ndof

        # remember reference nodal position
        self.X = X

        # remember material and geometric properties
        self.params = params

        # create an all zeros self.U
        self.U = []
        for X in self.X:
            self.U.append(np.zeros_like(X))

        # call additional element specific initialization steps
        self.init()

        # compute force and stiffness
        self.compute()

    def init(self):
        if not 'E' in self.params.keys():
            self.params['E'] = 1.0

        if not 'nu' in self.params.keys():
            self.params['nu'] = 0.0

    def setDisp(self, U):
        self.U = U
        self.compute()

    def getFe(self):
        return self.force

    def getKe(self):
        return self.stiffness

    def getFeAsMatrix(self):
        # return the internal force vector as nx1 matrix
        return np.matrix(np.hstack(self.force)).T

    def getKeAsMatrix(self):
        # return the stiffness matrix as nxn matrix
        return np.matrix(np.vstack([np.hstack(row) for row in self.stiffness]))

    def compute(self):

        # compute strain

        # compute stress

        # compute internal force vector
        self.force = copy.deepcopy(self.U)

        # compute tangent stiffness matrix
        self.stiffness = []
        for ndI in self.X:
            row = []
            for ndJ in self.X:
                KIJ = np.zeros((len(ndI),len(ndJ)))
                row.append(KIJ)
            self.stiffness.append(row)


# defining main execution procedure

def main():
    X1 = np.array([0.,0.])
    X2 = np.array([1.,0.])
    e = Element((X1,X2))
    print('U = ',e.U)
    print('---')
    print('Fe = ',e.getFe())
    print('Kte = \n',e.getKe())
    print('---')
    print('Fe = ',e.getFeAsMatrix())
    print('Kte = \n',e.getKeAsMatrix())


# main execution ****************************************

if __name__ == "__main__":
    main()
    sys.exit(0)
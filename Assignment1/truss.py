"""
class implementation of a 2D/3D truss class

creator: Peter Mackenzie-Helnwein
date: Oct 29, 2019
revised: April 5, 2020
"""

# import libraries
import numpy as np

# the element master class
from Element import *


# class definition
class TrussElement(Element):
    """
    !class: TrussElement

    variables:
        self.N0  ... unit normal vector from 1 to 2
        self.n  ... unit normal vector from 1 to 2

    inherited variables:
        self.nnode ......... number of nodes per element
        self.ndof .......... number of degrees of freedom per node
        self.X = (X1,X2) ... tuple of nodal position vectors (as np.array)
        self.U ............. array of nodal displacement vectors
        self.force ......... internal forces at nodes
        self.stiffness ..... the stiffness matrix

    overloaded methods:
        def init(self) ....... element specific initialization steps
        compute(self) ...... does the actual computation

    inherited methods:
        __init__(self, X, params)
        setDisp(self, U) ... U is an array of nodal displacement vectors
        getFe(self) ........ return the internal force vector as array of nodal vectors
        getKe(self) ........ return the stiffness matrix as array of nodal matrices
        getFeAsMatrix(self) .. return the internal force vector as nx1 matrix
        getKeAsMatrix(self) .. return the stiffness matrix as nxn matrix
    """

    def init(self):

        # make sure all needed parameters exist
        if not 'E' in self.params.keys():
            self.params['E'] = 1.0

        if not 'A' in self.params.keys():
            self.params['A'] = 1.0

        # compute reference base vectors
        L0vec = self.X[1] - self.X[0]
        L02 = L0vec @ L0vec
        self.L0 = np.sqrt(L02)
        self.N0 = L0vec / self.L0


    def compute(self):

        # compute deformed base vectors
        x1 = self.X[0] + self.U[0]
        x2 = self.X[1] + self.U[1]

        lvec = x2 - x1
        l2 = lvec @ lvec
        l  = np.sqrt(l2)
        n = lvec / l

        # compute strain
        stretch = l / self.L0
        strain = np.log(stretch)

        # compute internal force
        k = self.params['E'] * self.params['A']
        fe = k * strain

        # compute nodal force
        self.force = (-fe * n, fe * n)

        # compute tangent stiffness
        ke = (k - fe)/l * np.outer(n, n) + fe/l * np.identity(self.ndof)
        self.stiffness = np.array([[ke,-ke],[-ke,ke]])


# defining main execution procedure

def main():
    # create a demo element
    X1 = np.array([0.,0.])
    X2 = np.array([2.,1.])
    e = TrussElement((X1,X2))

    # undeformed state
    print('U = ',e.U)
    print('---')
    print('Fe = ',e.getFe())
    print('Kte = \n',e.getKe())
    print('---')
    print('Fe = ',e.getFeAsMatrix())
    print('Kte = \n',e.getKeAsMatrix())

    # now set some displacement and check out changes
    e.setDisp((np.array([0.0,0.0]),np.array([0.1,0.0])))
    print('=======')
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

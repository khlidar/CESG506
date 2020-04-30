"""
class implementation of a 2D rotation class

creator: Peter Mackenzie-Helnwein
date: Oct 29, 2019
revised: April 5, 2020
"""

# import libraries
import sys
import numpy as np

# the element master class
from Element import *


# class definition
class PlateElement(Element):
    """
    !class: PlateElement

    this is still the linear element

    variables:
        self.C  ... material stiffness tensor
        self.A  ... element area
        self.g1 ... covariant base vector 1
        self.g2 ... covariant base vector 2
        self.B1 ... kinematic matrix for node 1
        self.B2 ... kinematic matrix for node 2
        self.B3 ... kinematic matrix for node 3

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

        X1 = self.X[0]
        X2 = self.X[1]
        X3 = self.X[2]

        self.stress = np.zeros(3)

        if not 'E' in self.params.keys():
            self.params['E'] = 1.0

        if not 'nu' in self.params.keys():
            self.params['nu'] = 0.0

        E  = self.params['E']
        nu = self.params['nu']

        self.C = np.array([[1.,nu,0.],
                           [nu,1.,0.],
                           [0.,0.,(1-nu)/2.]]) * E/(1.-nu*nu)

        # compute base vectors
        self.g1 = X1 - X3
        self.g2 = X2 - X3

        # compute metric
        G = np.array([[np.dot(self.g1,self.g1),np.dot(self.g1,self.g2)],
                      [np.dot(self.g2,self.g1),np.dot(self.g2,self.g2)]])

        # compute jacobian and area
        J = np.sqrt( np.linalg.det(G) )
        self.A = J/2.

        # compute dual base vectors
        H = np.linalg.inv(G)
        h1 = H[0][0] * self.g1 + H[0][1] * self.g2
        h2 = H[1][0] * self.g1 + H[1][1] * self.g2

        self.B1 = np.array([[h1[0],0.0],[0.0,h1[1]],[h1[1],h1[0]]])
        self.B2 = np.array([[h2[0],0.0],[0.0,h2[1]],[h2[1],h2[0]]])
        self.B3 = -(self.B1 + self.B2)

    def compute(self):

        U1 = self.U[0]
        U2 = self.U[1]
        U3 = self.U[2]

        strain =  self.B1 @ (U1 - U3)
        strain += self.B2 @ (U2 - U3)
        stress = self.C @ strain


        self.force = (
             stress @ self.B1 * self.A,
             stress @ self.B2 * self.A,
             stress @ self.B3 * self.A
        )

        K11 = self.B1.T @ self.C @ self.B1 * self.A
        K12 = self.B1.T @ self.C @ self.B2 * self.A
        K13 = self.B1.T @ self.C @ self.B3 * self.A
        K22 = self.B2.T @ self.C @ self.B2 * self.A
        K23 = self.B2.T @ self.C @ self.B3 * self.A
        K33 = self.B3.T @ self.C @ self.B3 * self.A

        self.stiffness = ( (K11  ,K12  ,K13),
                           (K12.T,K22  ,K23),
                           (K13.T,K23.T,K33) )



# defining main execution procedure

def main():
    # create a demo element
    X1 = np.array([0.,0.])
    X2 = np.array([1.,0.])
    X3 = np.array([0.,1.])
    e = PlateElement((X1,X2,X3))

    # undeformed state
    print('U = ',e.U)
    print('---')
    print('Fe = ',e.getFe())
    print('Kte = \n',e.getKe())
    print('---')
    print('Fe = ',e.getFeAsMatrix())
    print('Kte = \n',e.getKeAsMatrix())

    # now set some displacement and check out changes
    e.setDisp((np.array([0.0,0.0]),np.array([0.1,0.0]),np.array([0.0,0.0])))
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


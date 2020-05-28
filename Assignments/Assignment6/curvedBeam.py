"""
class implementation of a 2D/3D curved beam

creator: Kristinn Hlíðar Grétarsson
date: 2020-05-15

"""

# import libraries
import numpy as np

# the element master class
from Element import *
from GaussianQuadrature import *
from numpy import array, vstack, append, hstack, outer, ix_, zeros
from scipy import integrate
import time



# class definition
class CurvedBeam(Element):
    """
    !class: CurvedBeam

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

        if not 'I' in self.params.keys():
            self.params['I'] = 1.0

        # compute reference base vectors
        L0vec = self.X[1] - self.X[0]
        L02 = L0vec @ L0vec
        self.L0 = np.sqrt(L02)
        self.N0 = L0vec / self.L0

        # Compute horizontal length
        x1 = self.X[0]
        x2 = self.X[1]
        self.lx = x2[0] - x1[0]

        self.U = array([zeros(3), zeros(3)])

    def update(self):
        '''
        Author: Kristinn Hlíðar
        Date: 2020.04.19
        Description: Created to rerun length vector for when start end end locations are updated
        :return:
        '''
        # compute reference base vectors
        L0vec = self.X[1] - self.X[0]
        L02 = L0vec @ L0vec
        self.L0 = np.sqrt(L02)
        self.N0 = L0vec / self.L0

        # Compute horizontal length
        x1 = self.X[0]
        x2 = self.X[1]
        self.lx = x2[0] - x1[0]

        self.U = array([zeros(3), zeros(3)])

    def u_0(self, x):
        u1 = self.U[0]
        u2 = self.U[1]

        Nu, dNu = self.Nu(x)

        u = array([u1[0], u2[0]])

        u0 = Nu.dot(u)

        du0 = dNu.dot(u)

        return u0, du0

    def h_func(self, x):
        lx = self.lx

        Nu, dNu = self.Nu(x)

        xi = self.X[0]
        xj = self.X[1]
        y = array([xi[1], xj[1]])

        h = (1 - x/lx)*xi[1] + (x/lx)*xj[1]
        h = Nu.dot(y)
        dh = (-1/lx)*xi[1] + (1/lx)*xj[1]
        dh = dNu.dot(y)
        return h, dh

    def v_func(self, x):
        # Get vertical length and calculate Xi ratio (x/lx)
        lx = self.lx
        r = x/lx

        Nv, dNv, ddNv = self.Nv(x)

        # Get displacement of element
        u1 = self.U[0]
        u2 = self.U[1]
        u = array([u1[1], u1[2], u2[1], u2[2]])
        '''
        v1 = (1 - 3*r**2 + 2*r**3) * u1[1]
        v2 = lx*(r - 2*r**2 + r**3) * u1[2]
        v3 = (3*r**2 - 2*r**3) * u2[1]
        v4 = lx*(r**3 - r**2) * u2[2]

        v = v1 + v2 + v3 + v4'''

        v = Nv.dot(u)
        '''
        dv1 = (-3 * r ** 2*(2/x) + 2 * r ** 3*(3/x)) * u1[1]
        dv2 = lx * (r*(1/x) - 2 * r ** 2*(2/x) + r ** 3*(3/x)) * u1[2]
        dv3 = (3 * r ** 2*(2/x) - 2 * r ** 3*(3/x)) * u2[1]
        dv4 = lx * (r ** 3*(3/x) - r ** 2*(2/x)) * u2[2]

        dv = dv1 + dv2 + dv3 + dv4'''
        dv = dNv.dot(u)
        '''
        ddv1 = (-3 * r ** 2 * (2 / x**2) + 2 * r ** 3 * (6 / x**2)) * u1[1]
        ddv2 = lx * (- 2 * r ** 2 * (2 / x**2) + r ** 3 * (6 / x**2)) * u1[2]
        ddv3 = (3 * r ** 2 * (2 / x**2) - 2 * r ** 3 * (6 / x**2)) * u2[1]
        ddv4 = lx * (r ** 3 * (6 / x**2) - r ** 2 * (2 / x**2)) * u2[2]

        ddv = ddv1 + ddv2 + ddv3 + ddv4'''
        ddv = ddNv.dot(u)

        return v, dv, ddv

    def eps_0(self, x):
        v, dv, ddv = self.v_func(x)
        h, dh = self.h_func(x)
        u0, du0 = self.u_0(x)

        eps = du0 + dh*dv + 0.5 * dv**2
        return eps

    def Force(self, x):
        eMod = self.params['E']
        area = self.params['A']
        eps = self.eps_0(x)

        f = eMod*area*eps

        return f

    def Moment(self, x):
        eMod = self.params['E']
        I = self.params['I']
        v, dv, phi = self.v_func(x)

        moment = eMod*I*phi
        return moment

    def Nv(self, x):
        # Get vertical length and calculate Xi ratio (x/lx)
        lx = self.lx
        r = x / lx

        v1 = (1 - 3 * r ** 2 + 2 * r ** 3)
        v2 = lx * (r - 2 * r ** 2 + r ** 3)
        v3 = (3 * r ** 2 - 2 * r ** 3)
        v4 = lx * (r ** 3 - r ** 2)

        Nv = array([v1, v2, v3, v4])

        dv1 = (-3 * r ** 2 * (2 / x) + 2 * r ** 3 * (3 / x))
        dv2 = lx * (r * (1 / x) - 2 * r ** 2 * (2 / x) + r ** 3 * (3 / x))
        dv3 = (3 * r ** 2 * (2 / x) - 2 * r ** 3 * (3 / x))
        dv4 = lx * (r ** 3 * (3 / x) - r ** 2 * (2 / x))

        dNv = array([dv1, dv2, dv3, dv4])

        ddv1 = (-3 * r ** 2 * (2 / x ** 2) + 2 * r ** 3 * (6 / x ** 2))
        ddv2 = lx * (- 2 * r ** 2 * (2 / x ** 2) + r ** 3 * (6 / x ** 2))
        ddv3 = (3 * r ** 2 * (2 / x ** 2) - 2 * r ** 3 * (6 / x ** 2))
        ddv4 = lx * (r ** 3 * (6 / x ** 2) - r ** 2 * (2 / x ** 2))

        ddNv = array([ddv1, ddv2, ddv3, ddv4])

        return Nv, dNv, ddNv

    def Nu(self, x):
        lx = self.lx

        u0_1 = (1 - x / lx)
        u0_2 =(x / lx)
        Nu = array([u0_1, u0_2])

        du0_1 = (-1 / lx)
        du0_2 = (1 / lx)
        dNu = array([du0_1, du0_2])

        return Nu, dNu

    def residual_x(self, x):
        force = self.Force(x)
        moment = self.Moment(x)
        v, dv, ddv = self.v_func(x)
        h, dh = self.h_func(x)

        lx = self.lx
        r = x/lx

        Nu, dNu = self.Nu(x)
        Nv, dNv, ddNv = self.Nv(x)

        r1 = force*dNu
        r2 = moment*ddNv + force*(dh+dv)*dNv
        r = hstack((r1, r2))

        r1 = r[[0, 2, 3]]
        r2 = r[[1, 4, 5]]

        r = array([r1, r2])

        return r

    def Residual(self, gpoints=4):
        lx = self.lx
        #I = integrate.quad_vec(self.residual_x, 0, lx)
        points, weights = Gaussian(gpoints)
        residual = array([zeros(3), zeros(3)])
        for i in range(0, gpoints):
            xi = lx * (1 + points[i]) / 2
            residual = residual + self.residual_x(xi) * lx * 0.5 * weights[i]

        return residual

    def Stiffness_x(self, x):
        Nu, dNu = self.Nu(x)
        Nv, dNv, ddNv = self.Nv(x)

        eMod = self.params['E']
        A = self.params['A']
        I = self.params['I']

        force = self.Force(x)
        v, dv, ddv = self.v_func(x)
        h, dh = self.h_func(x)

        stiff11 = eMod*A*outer(dNu, dNu)
        stiff12 = eMod*A*(dh + dv) * outer(dNu, dNv)
        stiff21 = eMod*A*(dh + dv) * outer(dNv, dNu)
        stiff22 = (eMod*A*(dh + dv)**2 + force) * outer(dNv, dNv) \
                    + eMod*I*outer(ddNv, ddNv)

        stiffness = vstack((hstack((stiff11, stiff12)), hstack((stiff21, stiff22))))

        stiff11 = stiffness[np.ix_([0, 2, 3], [0, 2, 3])]
        stiff12 = stiffness[np.ix_([0, 2, 3], [1, 4, 5])]
        stiff21 = stiffness[np.ix_([1, 4, 5], [0, 2, 3])]
        stiff22 = stiffness[np.ix_([1, 4, 5], [1, 4, 5])]

        stiffness = array([[stiff11, stiff12], [stiff21, stiff22]])

        return stiffness

    def Stiffness(self, gpoints=4):
        #start_time = time.time()
        lx = self.lx
        #I = integrate.quad_vec(self.Stiffness_x, 0, lx)
        #print("--- %s seconds element ---" % (time.time() - start_time))

        #start_time = time.time()
        stiff = array([[zeros((3, 3)), zeros((3, 3))], [zeros((3, 3)), zeros((3, 3))]])
        points, weights = Gaussian(gpoints)
        for i in range(0, gpoints):
            xi = lx * (1 + points[i]) / 2
            stiff = stiff + self.Stiffness_x(xi) * lx * 0.5 * weights[i]
        #print("--- %s seconds My Gaussian ---" % (time.time() - start_time))
        return stiff


    def compute(self):
        #start_time = time.time()
        if self.lx:
            # compute nodal force
            self.force = self.Residual()

            # compute tangent stiffness
            self.stiffness = self.Stiffness()
            pass
        #print("--- %s seconds Compute ---" % (time.time() - start_time))

# defining main execution procedure

def main():
    # create a demo element
    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}
    X1 = np.array([0., 0.])
    X2 = np.array([2., 0.])
    e = CurvedBeam((X1, X2), para)

    # undeformed state
    print('U = ',e.U)
    print('---')
    print('Fe = ',e.getFe())
    print('Kte = \n',e.getKe())
    print('---')
    print('Fe = ',e.getFeAsMatrix())
    print('Kte = \n',e.getKeAsMatrix())

    # now set some displacement and check out changes
    e.setDisp((np.array([0.0, 0.0, 0.]), np.array([-0.0025015644, 0.1, 0.])))
    print('=======')
    print('U = ',e.U)
    print('---')
    print('Fe = ',e.getFe())
    print('Kte = \n',e.getKe())
    print('---')
    print('Fe = ',e.getFeAsMatrix())
    print('Kte = \n',e.getKeAsMatrix())
    print(f'u_0(1) = {e.u_0(1)}')
    print(f'u_func(0.5) = {e.h_func(0.5)}')
    print(f'v_func(0.5) = {e.v_func(0.5)}')
    print(f'h_func(0.5) = {e.h_func(0.5)}')
    print(f'eps(0.5) = {e.eps_0(0.5)}')
    print(f'Force(0.5) = {e.Force(0.5)}')
    print(f'residual_x(0.5) = {e.residual_x(0.5)}')
    print(f'Residual = {e.Residual()}')
    print(f'Stiffness_x = {e.Stiffness_x(0.5)}')
    print(f'Stiffness = \n {e.Stiffness()}')

# main execution ****************************************

if __name__ == "__main__":
    main()
    sys.exit(0)

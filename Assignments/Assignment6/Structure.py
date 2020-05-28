'''
Created by: Kristinn Hlíðar Grétarsson
Date:       2020.05.19

Description:
    FEM structure program that uses FEM elements from Peter to construct and calculate non-linear problems

Functions:

'''

# Imports
from Element import *
from truss import *
from numpy import array, append, vstack, zeros, linalg, dot, hstack, transpose, newaxis, sqrt, linspace, inf, zeros_like
import matplotlib.pyplot as plt
from curvedBeam import *
from curvedBeam_V2 import *
import time

# Functions and classes
class struct(object):
    def __init__(self):
        self.nodes = array(['no.', 'State', '[x, y, z]'])
        self.free = array([], int)
        self.support = array([], int)

        self.ele = array(['no.', 'i-node', 'j-node', 'element'])
        self.ktot = array([])
        self.rtot = array([])

        self.loads = array(['no.', 'node', '[x,y,z]'])
        self.P = array([])

        self.u = array([])

        self.l = linspace(0, 1, 6)

        self.ndof = 0 # Number of degrees of freedom
        self.space = 0 # 1D, 2D or 3D space



    def addNode(self, nu, state=array([0, 0, 0]), loc=array([0, 0, 0])):
        '''
        :param nu: number of node, to be changed out for automatic counter
        :param state: State of node 0: support, 1: free
        :param loc: Location of node numpy.array
        :return: [nothing]
        '''
        temp = array([nu, state, loc])
        self.nodes = vstack((self.nodes, temp))

        # set number of degrees of freedom for structure
        ndof = len(loc)
        ndof = len(state) # ATT potential solution to nDoF
        if not self.ndof:
            self.ndof = ndof
        elif not self.ndof == ndof:
            print('Warning NDoF for current node not same as for previous nodes')

        # set dimension structure
        dim = len(loc)
        if not self.space:
            self.space = dim
        elif not self.space == dim:
            print('Warning NDoF for current node not same as for previous nodes')

        # add to list of free/supported movement
        count = 0
        for i in state:
            temp = (nu - 1) * ndof + count
            temp = int(temp)
            if i:
                self.free = append(self.free, temp)
            else:
                self.support = append(self.support, temp)
            count += 1

        self.u = zeros(self.getNoNodes()*ndof)

    def addEle(self, nu, iNode, jNode, ele:Element):
        '''
            Adds element to the list of elements in structure
        :param nu: Number of element, to be changed out for automatic counter
        :param iNode: Starting node of element
        :param jNode: End node of element
        :param ele: Element
        :return:
        '''
        temp = array([nu, iNode, jNode, ele])
        self.ele = vstack((self.ele, temp))

    def addLoad(self, nu, node, load:array([])):
        temp = array([nu, node, load])
        self.loads = vstack((self.loads, temp))
        self.updateLoads()

    def calcK(self):
        '''
        Calculates the stiffness matrix of the structure
        Stores the new stiffness matrix in the self.ktot variable
        :return:
        '''
        nn = self.getNoNodes()
        n = self.ndof
        n = 3 # ATT find better solution
        ktot = zeros((nn*n, nn*n))
        for ele in self.ele:
            if not ele[0] == 'no.':
                inode = ele[1]
                jnode = ele[2]

                i = (inode - 1) * n
                j = (jnode - 1) * n

                loci = self.nodes[inode, 2]
                locj = self.nodes[jnode, 2]

                kele = ele[3].getKe()

                ktot[i:i + n, i:i + n] += kele[0, 0]
                ktot[i:i + n, j:j + n] += kele[0, 1]
                ktot[j:j + n, i:i + n] += kele[1, 0]
                ktot[j:j + n, j:j + n] += kele[1, 1]

        self.ktot = ktot

    def calcR(self):
        '''
        Calculates the current forces due to deformation of the structure
        Stores the forces in the variable self.rtot
        :return:
        '''
        nn = self.getNoNodes()
        n = self.ndof
        n = 3 # ATT Find better solution
        Rtot = zeros((nn * n))
        for ele in self.ele:
            if not ele[0] == 'no.':
                inode = ele[1]
                jnode = ele[2]

                i = (inode - 1) * n
                j = (jnode - 1) * n

                loci = self.nodes[inode, 2]
                locj = self.nodes[jnode, 2]

                Rele = ele[3].getFe()

                Rtot[i:i + n] += Rele[0]
                Rtot[j:j + n] += Rele[1]

        self.rtot = Rtot

        return Rtot

    def partitionK(self):
        '''
        Partitions the stiffness matrix according to DoF
        '''
        self.calcK()

        self.kff = self.ktot[self.free[:, None], self.free]
        self.kfs = self.ktot[self.free[:, None], self.support]
        self.ksf = self.ktot[self.support[:, None], self.free]
        self.kss = self.ktot[self.support[:, None], self.support]

    def updateElements(self):
        for ele in self.ele:
            if not ele[0] == 'no.':
                inode = ele[1]
                jnode = ele[2]

                loci = self.nodes[inode, 2]
                locj = self.nodes[jnode, 2]

                ele[3].setLoc((loci, locj))

    def updateLoads(self):
        n = self.getNoNodes()
        ndof = self.ndof
        # Create load array with as many DoF the structure has
        self.P = zeros(n*ndof)

        # Loop through all loads on structure and add to total load array
        for load in self.loads:
            if not load[0] == 'no.':
                node = load[1]
                for i in range(ndof):
                    self.P[(node-1)*ndof + i] += load[2][i]
                    pass

    def getNodes(self, show=False):
        if show:
            print(self.nodes)
        return self.nodes

    def getNoNodes(self):
        output = 0
        for nodes in self.nodes:
            output += 1
        return output -1

    def getElements(self, show=False):
        if show:
            print(self.ele)
        return self.ele

    def getElementLoc(self):
        output = array(['no.', 'i-node', 'j-node'])
        for ele in self.ele:
            if not ele[0] == 'no.':
                inode = self.nodes[ele[1], 2]
                jnode = self.nodes[ele[2], 2]
                temp = array([ele[0], inode, jnode])
                output = vstack((output, temp))
        return output

    def solve(self):
        tol = 1e-12

        self.resetDisp()

        fdof = len(self.free)
        u = zeros(fdof)
        du = zeros(fdof)

        self.calcR()

        for lpf in self.l:
            rf = self.rtot[self.free]

            pf = self.P[self.free] * lpf

            error = 0
            for i in range(len(rf)):
                error += abs(rf[i] - pf[i])

            count = 0
            while error > tol:
                self.partitionK()
                kff = self.kff

                du = linalg.solve(kff, pf - rf)

                u += du

                self.u[self.free] = u

                self.displaceEle()

                rf = self.rtot[self.free]
                error = 0
                for i in range(len(rf)):
                    error += abs(rf[i] - pf[i])**2
                error = sqrt(error)

                count += 1
                if count > 20:
                    print('Failed to converge to solution')
                    print(f'Load step {lpf}')
                    break

        return self.u

    def solveDis(self, disControl:array([]), g=1, reset=True, la=1):
        tol = 1e-10

        # initialize load factor
        #la = 1

        # make sure to reset displacements
        if reset:
            self.resetDisp()

        # Store initial displacement
        u_ini = self.u

        # Get which node user wants to displace and in what direction
        disNode = disControl[0]
        load_dir = disControl[1]

        # Create load and displacement vector
        pref = zeros(self.getNoNodes() * self.ndof)
        ek = zeros(self.getNoNodes() * self.ndof)
        du = self.u[self.free]

        # Add reference force and displacement direction to correct location in
        i = (disNode-1) * self.ndof
        pref[i:i+self.ndof] = load_dir
        pref = pref[self.free]

        ek[i:i + self.ndof] = load_dir
        ek = ek[self.free]
        ek = abs(ek)

        # Calculate internal forces
        self.calcR()
        rf = self.rtot[self.free]

        # Calculate force
        pf = pref * la

        hasfailed = False

        for lpf in self.l:

            error = 1
            count = 0
            while error > tol:
                # partition and get stiffness matrix
                self.partitionK()
                kff = self.kff

                du1 = linalg.solve(kff, pref)
                du0 = linalg.solve(kff, pf - rf)

                #test, test1 = linalg.solve(kff, (pref, pf-rf))

                dla = -(dot(ek, du) - g * lpf + dot(ek, du0))/(dot(ek, du1))
                la += dla

                # Displacement
                du += du0 + dla*du1

                # Set displacement
                self.u[self.free] = du
                self.displaceEle()

                # Calculate forces
                rf = self.calcR()[self.free]

                # Calculate force in direction of displacement
                pf = pref * la

                error = 0
                for i in range(len(rf)):
                    error += abs(rf[i] - pf[i])**2
                error = sqrt(error)

                count += 1

                if count > 20:
                    print(f'Failed to converge dis. control. Dis. = {g}')
                    print(f'det(kff) = {linalg.det(kff)}')
                    #self.u = u_ini
                    #du = u_ini[self.free]
                    #self.displaceEle()
                    #la = 1
                    hasfailed = True

                    break

        return la, self.u, hasfailed

    def solveArc(self, s_lim=inf, load=inf, nodeDis=array([0]), reset=True, la=0.1, alfa=0, data=False):
        '''

        :param disNode: Which node is being used for control
        :param s_lim: Desired total archlength
        :param reset: Should structure displacement be reset before solving
        :param la: Initial load factor used for first step
        :param alfa:
        :return:
        '''

        # Equilibrium tolarance
        tol = 1e-8

        # make sure to reset displacements
        if reset:
            self.resetDisp()

        # Set last displacement and load factor
        u_n = self.u[self.free]
        la_n = 0 # Come back to this, user should be able to set first KHG

        # Create load and displacement vector
        pref = self.P[self.free]
        du = self.u[self.free]
        u_hist = self.u
        la_hist = array([la_n])
        s_hist = 0
        R_hist = 0
        det_hist = 0

        # Calculate internal forces
        self.calcR()
        rf = self.rtot[self.free]

        hasfailed = False

        # Total arc length and delta_s initialized
        s_tot = 0
        ds = 0

        # node displacement control setup
        if nodeDis[0]:
            # node with target displacement
            node = nodeDis[0]

            # Vector to locate target deflection
            ek = zeros(self.getNoNodes() * self.ndof)


            # Add reference force and displacement direction to correct location in
            i = (node - 1) * self.ndof
            print(nodeDis[1])
            ek[i:i + self.ndof] = nodeDis[1]
            ek = ek[self.free]

            for t in nodeDis[1]:
                if t:
                    nodeTarget = t

            ek = ek / linalg.norm(ek)
            ek = abs(ek)


            pass

        else:
            ek = zeros_like(du)
            nodeTarget = -inf


        while s_tot < s_lim and la <= load and nodeTarget - dot(ek, du) < 0:

            error = 1
            count = 0

            print('-+-+-+-+-+-+-+-+-+-+-')

            while error > tol:
                #start_time = time.time()
                # partition and get stiffness matrix
                self.partitionK()
                kff = self.kff
                #print("--- %s seconds update [K] ---" % (time.time() - start_time))

                # if delta_s hasn't been calculated then this is the first iteration and we calculate delta_s
                if count == 0 and not ds:
                    du = linalg.solve(kff, pref*la)
                    ds2 = dot(du, du) + alfa * la**2 * dot(pref, pref)
                    ds = sqrt(ds2)

                    # Set displacement
                    self.u[self.free] = du
                    self.displaceEle()

                    # Calculate forces
                    rf = self.calcR()[self.free]

                    #R0 = rf - pref * la
                    #R0 = append(R0, 0)
                    #R0 = linalg.norm(R0)
                elif count == 0:
                    du = 2*u_n - u_n1
                    la = 2*la_n - la_n1

                    #g = dot(du - u_n, du - u_n) + alfa * (la - la_n) ** 2 * dot(pref, pref) - ds2
                    #R0 = rf - pref * la
                    #R0 = append(R0, g)
                    #R0 = linalg.norm(R0)
                    pass

                g = dot(du-u_n, du-u_n) + alfa*(la - la_n)**2 * dot(pref, pref) - ds2

                #start_time = time.time()
                du1 = linalg.solve(kff, pref)
                du0 = linalg.solve(kff, pref*la - rf)
                residual = pref*la - rf
                #print("--- %s seconds solve displacement ---" % (time.time() - start_time))

                dla = -(g + dot(2*(du-u_n), du0)) / (dot(2*(du-u_n), du1) + 2*alfa*(la-la_n)*dot(pref, pref))
                la += dla

                # Displacement
                du += du0 + dla * du1

                #start_time = time.time()
                # Set displacement
                self.u[self.free] = du
                self.displaceEle()
                #print("--- %s seconds set disp. ---" % (time.time() - start_time))

                #start_time = time.time()
                # Calculate forces
                rf = self.calcR()[self.free]
                residual = pref * la - rf
                #print("--- %s seconds calc forces ---" % (time.time() - start_time))

                #start_time = time.time()
                err = rf - pref*la
                err = append(err, g)
                error = linalg.norm(err)

                # testing normalised error
                #error = min(error / R0, error)
                print(error)
                #print("--- %s seconds calc error ---" % (time.time() - start_time))

                # Before moving on to next step set last convergence
                if error < tol:
                    u_n1 = u_n
                    la_n1 = la_n

                    u_n = du
                    la_n = la

                    s_tot += ds

                    u_hist = vstack((u_hist, self.u))
                    la_hist = vstack((la_hist, la))

                    if data:
                        eigv = linalg.eigvals(self.kff)
                        #det_hist = vstack((det_hist, linalg.det(self.kff)))
                        det_hist = vstack((det_hist, min(eigv)))
                        pass

                if data:
                    s_hist = vstack((s_hist, s_tot))
                    R_hist = vstack((R_hist, error))


                count += 1
                if count > 20:
                    print(f'Failed to converge dis. control. Dis. = {s_tot}')
                    print(f'det(kff) = {linalg.det(kff)}')
                    # self.u = u_ini
                    # du = u_ini[self.free]
                    # self.displaceEle()
                    # la = 1
                    hasfailed = True

                    break


        if data:
            return u_hist, la_hist, s_hist, R_hist, det_hist
        else:
            return u_hist, la_hist

        pass

    def displacementpath(self, disControl:array([]), g=[0, 1], steps=20):
        dis = linspace(g[0], g[1], steps)
        upath = zeros(len(self.free))
        load = [0]
        load1 = 1

        for d in dis:
            load1, u, failed = self.solveDis(disControl, d, reset=False, la=load1)
            if not failed:
                upath = vstack((upath, u[self.free]))
                load.append(load1)
            else:
                load1 = load[len(load)-1]

        return load, upath

    def plotStructure(self, deformed=False, mark=False):
        if self.space == 2:
            fig = plt.figure()
            ax = fig.add_subplot(111)

            for node in self.nodes:
                if not node[0] == 'no.':
                    ax.scatter(node[2][0], node[2][1], color='b', marker='o')

                    if mark:
                        plt.text(node[2][0], node[2][1], node[0], color='b')

            shape = self.ele.shape
            if not len(shape) == 1:
                for ele in self.ele:
                    if not ele[0] == 'no.':
                        i = ele[1]
                        j = ele[2]

                        x = [self.nodes[i][2][0], self.nodes[j][2][0]]
                        y = [self.nodes[i][2][1], self.nodes[j][2][1]]

                        ax.plot(x, y, 'b')

                        if mark:
                            plt.text((x[0]+x[1])/2, (y[0]+y[1])/2, ele[0], color='blue',
                                     bbox=dict(facecolor='none', edgecolor='blue'))

                if deformed:
                    for node in self.nodes:
                        if not node[0] == 'no.':
                            loc = (node[0] - 1) * self.ndof
                            x = node[2][0] + self.u[loc]
                            y = node[2][1] + self.u[loc + 1]
                            ax.scatter(x, y, color='r', marker='o')

                    for ele in self.ele:
                        if not ele[0] == 'no.':
                            i = ele[1]
                            j = ele[2]

                            loc = (i - 1) * self.ndof
                            xi = self.nodes[i][2][0] + self.u[loc]
                            yi = self.nodes[i][2][1] + self.u[loc + 1]

                            loc = (j - 1) * self.ndof
                            xj = self.nodes[j][2][0] + self.u[loc]
                            yj = self.nodes[j][2][1] + self.u[loc + 1]

                            x = [xi, xj]
                            y = [yi, yj]

                            ax.plot(x, y, 'r')
            plt.show()
        elif self.space == 3:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            for node in self.nodes:
                if not node[0] == 'no.':
                    ax.scatter(node[2][0], node[2][1], node[2][2], color='b', marker='o')

            for ele in self.ele:
                if not ele[0] == 'no.':
                    i = ele[1]
                    j = ele[2]

                    x = [self.nodes[i][2][0], self.nodes[j][2][0]]
                    y = [self.nodes[i][2][1], self.nodes[j][2][1]]
                    z = [self.nodes[i][2][2], self.nodes[j][2][2]]

                    ax.plot(x, y, z, 'b')

            if deformed:
                for node in self.nodes:
                    if not node[0] == 'no.':
                        loc = (node[0] - 1) * self.ndof
                        x = node[2][0] + self.u[loc]
                        y = node[2][1] + self.u[loc + 1]
                        z = node[2][2] + self.u[loc + 2]

                        ax.scatter(x, y, z, color='r', marker='o')

                for ele in self.ele:
                    if not ele[0] == 'no.':
                        i = ele[1]
                        j = ele[2]

                        loc = (i - 1) * self.ndof
                        xi = self.nodes[i][2][0] + self.u[loc]
                        yi = self.nodes[i][2][1] + self.u[loc + 1]
                        zi = self.nodes[i][2][2] + self.u[loc + 2]

                        loc = (j - 1) * self.ndof
                        xj = self.nodes[j][2][0] + self.u[loc]
                        yj = self.nodes[j][2][1] + self.u[loc + 1]
                        zj = self.nodes[j][2][2] + self.u[loc + 2]

                        x = [xi, xj]
                        y = [yi, yj]
                        z = [zi, zj]

                        ax.plot(x, y, z, 'r')

            plt.show()

    def displaceEle(self, nodes=[], du=array([])):
        ndof = self.ndof
        n = self.getNoNodes()

        if len(du):
            self.u[nodes] = du

        for ele in self.ele:

            if not ele[0] == 'no.':

                inode = ele[1]
                jnode = ele[2]

                i = (inode-1)*ndof
                j = (jnode-1)*ndof

                ui = self.u[i:i+ndof]
                uj = self.u[j:j+ndof]
                #start_time = time.time()
                ele[3].setDisp((ui, uj))
                #print("--- %s seconds disp. displacement ---" % (time.time() - start_time))
        self.calcR()

    def resetDisp(self):
        for ele in self.ele:
            if not ele[0] == 'no.':
                ndof = self.ndof
                ndof = 3 # ATT Find better solution

                inode = ele[1]
                jnode = ele[2]

                loci = self.nodes[inode, 2]
                locj = self.nodes[jnode, 2]

                ele[3].setDisp((zeros(ndof), zeros(ndof)))






# defining main execution procedure
def main():
    p2 = struct()

    nrNodes = 5
    length = 500
    h0 = length/10
    x = linspace(0, length, nrNodes)
    y = (4*h0*(x/length))*(1-(x/length))  # curvature of beam
    #y = zeros_like(x)
    for i in range(0, nrNodes):
        if i == 0:
            p2.addNode(i+1, array([0, 0, 1]), array([x[i], y[i]]))
        elif i == nrNodes-1:
            p2.addNode(i + 1, array([0, 0, 1]), array([x[i], y[i]]))
        else:
            p2.addNode(i+1, array([1, 1, 1]), array([x[i], y[i]]))


    para = {'E': 1.0,
            'A': 10000.0,
            'I': 100000.0,
            'nu': 0.0}

    for i in range(0, nrNodes - 1):
        p2.addEle(i+1, i+1, i+2, CurvedBeam_V2(params=para))



    p2.updateElements()
    #p2.getElements(True)




    #p2.calcK()
    #p2.calcR()

    p2.partitionK()

    p2.addLoad(1, 1, array([0., 0., 0.]))
    p2.addLoad(1, 3, array([0., -1., 0.]))
    start_time = time.time()
    upath, load = p2.solveArc(load=1, la=0.2, alfa=0.05)
    print("--- %s seconds ---" % (time.time() - start_time))
    print(upath)
    print(load)
    p2.plotStructure(deformed=True)

    start_time = time.time()
    p2.partitionK()
    print(p2.kff)
    print(p2.rtot)
    print("--- %s seconds ---" % (time.time() - start_time))

    '''
    
    

    upath, load = p2.solveArc(load=2)

    p2.plotStructure(True)

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    #ax1.scatter(upath[:, 0] + 5.5, upath[:, 1] + 3.75, upath[:, 2] + 0.5, color='b', label='$node 5$')
    #ax1.scatter(upath[:, 3] + 5.5, upath[:, 4] + 2.5, upath[:, 5] + 0.5, color='r', label='$node 6$')

    ax1.scatter(upath[:, 0], upath[:, 1], upath[:, 2], color='b', label='$node 5$')
    ax1.scatter(upath[:, 3], upath[:, 4], upath[:, 5], color='r', label='$node 6$')

    ax1.legend()

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(upath[:, 0], load, color='b', label='$u_x$')
    ax1.scatter(upath[:, 1], load, color='r', label='$u_y$')
    ax1.scatter(upath[:, 2], load, color='g', label='$u_z$')
    ax1.legend()

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(111)
    ax1.scatter(upath[:, 3], load, color='b', label='$u_x$')
    ax1.scatter(upath[:, 4], load, color='r', label='$u_y$')
    ax1.scatter(upath[:, 5], load, color='g', label='$u_z$')
    ax1.legend()

    plt.show()'''


    '''
    p2.addNode(1, array([0, 0]), array([0.0, 0.0]))
    p2.addNode(2, array([1, 1]), array([5.5, 0.5]))
    p2.addNode(3, array([0, 0]), array([9.5, 0.0]))

    para = {'E': 2100.0,
            'A': 1.0,
            'nu': 0.0}

    p2.addEle(1, 1, 2, TrussElement(params=para))
    p2.addEle(2, 2, 3, TrussElement(params=para))

    p2.updateElements() 

    p = p2.solveDisV2(array([2, array([0, 1])]), d=0.1)

    print(p)

    p2.addLoad(1, 2, array([0., p]))

    u = p2.solve()
    print(u)

    print(p2.rtot)'''

# main execution ****************************************
if __name__ == "__main__":
    main()
    sys.exit(0)
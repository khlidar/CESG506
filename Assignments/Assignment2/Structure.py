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
from numpy import array, append, vstack, zeros, linalg, dot, hstack, transpose, newaxis

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

        self.l = [0, 0.25, 0.5, 0.75, 0.9, 1]

        self.ndof = 0


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
        if not self.ndof:
            self.ndof = ndof
        elif not self.ndof == ndof:
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
        nn = self.getNoNodes()
        n = self.ndof
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
        nn = self.getNoNodes()
        n = self.ndof
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
        self.calcK()

        self.kff = self.ktot[self.free[:,None], self.free]
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
        self.P = zeros(n*ndof)

        for load in self.loads:
            if not load[0] == 'no.':
                node = load[1]
                for i in range(ndof):
                    self.P[(node-1)*ndof + i] = load[2][i]
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
        for ele in self.ele:
            if not ele[0] == 'no.':
                loc = ele[3].getLoc()

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
                    error += abs(rf[i] - pf[i])

                count += 1
                if count > 20:
                    print('Failed to converge to solution')
                    print(f'Load step {lpf}')
                    break

        return self.u

    def solveDisV2(self, disControl:array([]), d=1):
        tol = 1e-12

        self.resetDisp()

        self.resetDisp()
        disNode = disControl[0]
        load_dir = disControl[1]

        ptot = zeros(self.getNoNodes() * self.ndof)
        ek = zeros(self.getNoNodes() * self.ndof)

        i = (disNode - 1) * self.ndof
        ptot[i:i + self.ndof] = load_dir
        ek[i:i + self.ndof] = load_dir
        ek = ek[self.free]

        fdof = len(self.free)
        l = 1
        u = zeros(fdof)
        du = zeros(fdof)

        self.calcR()

        for lpf in [1]:
            rf = self.rtot[self.free]

            pf = ptot[self.free] * lpf

            error = 0
            for i in range(len(rf)):
                error += abs(rf[i] - pf[i])

            count = 0
            while error > tol:
                self.partitionK()
                kff = self.kff
                pft = transpose(pf[newaxis])
                kff = hstack((kff, pft))
                kff = vstack((kff, hstack((ek, array([0.])))))


                R = pf - rf
                g = d - dot(ek, u)

                R = append(R, g)

                du = linalg.solve(kff, R)

                l -= du[len(du)-1]
                u += du[0:len(du)-1]

                pf = ptot[self.free]*l

                self.u[self.free] = u

                self.displaceEle()

                rf = self.rtot[self.free]
                error = 0
                for i in range(len(rf)):
                    error += abs(rf[i] - pf[i])

                count += 1
                if count > 20:
                    print('Failed to converge to solution')
                    print(f'Load step {lpf}')
                    break

        return l

        pass

    def solveDis(self, disControl:array([]), g=1):
        tol = 1e-6

        self.resetDisp()
        disNode = disControl[0]
        load_dir = disControl[1]

        ptot = zeros(self.getNoNodes() * self.ndof)
        utot = zeros(self.getNoNodes() * self.ndof)

        i = (disNode-1) * self.ndof
        ptot[i:i+self.ndof] = load_dir
        utot[i:i + self.ndof] = load_dir

        pf = ptot[self.free]
        uf = utot[self.free]

        self.calcR()
        rf = self.rtot[self.free]

        self.partitionK()

        kff = self.kff

        u1 = linalg.solve(kff, pf)

        du0 = linalg.solve(kff, rf)

        la = (g + dot(uf, du0))/(dot(uf, u1))

        du = du0 + la*u1

        self.u[self.free] = du

        self.displaceEle()

        rtot = self.calcR()
        rf = rtot[self.free]
        pf = ptot[self.free]*la

        error = 0
        for i in range(len(rf)):
            error += abs(rf[i] - pf[i])

        while error > tol:
            self.partitionK()

            kff = self.kff

            u1 = linalg.solve(kff, pf)

            du0 = linalg.solve(kff, pf - rf)

            la = -(-g + dot(uf, du0)) / (dot(uf, u1))

            du = du0 + la * u1

            self.u[self.free] = du

            self.displaceEle()

            rtot = self.calcR()
            rf = rtot[self.free]
            #pf = ptot[self.free] * la
            pf = pf * la
            #pf = ptot[self.free]

            error = 0
            for i in range(len(rf)):
                error += abs(rf[i] - pf[i])

            pass




        print(rf)
        print(pf)


        return la

    def displaceEle(self):
        ndof = self.ndof
        n = self.getNoNodes()

        for ele in self.ele:
            if not ele[0] == 'no.':
                inode = ele[1]
                jnode = ele[2]

                i = (inode-1)*ndof
                j = (jnode-1)*ndof

                ui = self.u[i:i+ndof]
                uj = self.u[j:j+ndof]

                ele[3].setDisp((ui, uj))

        self.calcR()

    def resetDisp(self):
        for ele in self.ele:
            if not ele[0] == 'no.':
                inode = ele[1]
                jnode = ele[2]

                loci = self.nodes[inode, 2]
                locj = self.nodes[jnode, 2]

                ele[3].setDisp((zeros(self.ndof), zeros(self.ndof)))






# defining main execution procedure
def main():
    p2 = struct()


    p2.addNode(1, array([0, 0, 0]), array([0.0, 5.  , 0.]))
    p2.addNode(2, array([0, 0, 0]), array([9.5, 5.  , 0.]))
    p2.addNode(3, array([0, 0, 0]), array([0.0, 0.  , 0.]))
    p2.addNode(4, array([0, 0, 0]), array([9.5, 0.  , 0.]))
    p2.addNode(5, array([1, 1, 1]), array([5.5, 3.75, 0.5]))
    p2.addNode(6, array([1, 1, 1]), array([5.5, 2.50, 0.5]))


    para = {'E': 2100.0,
            'A': 1.0,
            'nu':0.0}

    p2.addEle(1, 1, 5, TrussElement(params=para))
    p2.addEle(2, 1, 6, TrussElement(params=para))
    p2.addEle(3, 2, 5, TrussElement(params=para))
    p2.addEle(4, 3, 6, TrussElement(params=para))
    p2.addEle(5, 4, 6, TrussElement(params=para))
    p2.addEle(6, 4, 5, TrussElement(params=para))
    p2.addEle(7, 5, 6, TrussElement(params=para))

    p2.updateElements()
    p2.getElements(True)

    p2.calcK()
    p2.calcR()
    p2.partitionK()

    p = p2.solveDisV2(array([5, array([0, 0, 1])]), d=0.01)

    print(p)

    p2.addLoad(1, 5, array([0., 0., p]))

    u = p2.solve()
    print(u)


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
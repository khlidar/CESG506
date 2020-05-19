
from Structure import *


def load(x):
    w = [-.99, -1.01]

    wx = (w[0] + (w[1] - w[0])/500*x)*0.01

    return wx

def Nv(x, lx):
    # Get vertical length and calculate Xi ratio (x/lx)
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


p5 = struct()

nrNodes = 65  # Number of nodes
length = 500  # Length between
h0 = length/10  # height at mid span
x = linspace(0, length, nrNodes)
y = (4*h0*(x/length))*(1-(x/length))  # curvature of beam
#y = zeros_like(x)
for i in range(0, nrNodes):
    if i == 0:
        p5.addNode(i + 1, array([0, 0, 1]), array([x[i], y[i]]))
    elif i == nrNodes-1:
        p5.addNode(i + 1, array([0, 0, 1]), array([x[i], y[i]]))
    else:
        p5.addNode(i + 1, array([1, 1, 1]), array([x[i], y[i]]))

para = {'E': 1.0,
        'A': 10000.0,
        'I': 100000.0,
        'nu': 0.0}

for i in range(0, nrNodes - 1):
    p5.addEle(i + 1, i + 1, i + 2, CurvedBeam(params=para))

# Assign loads
lx = length / (nrNodes - 1)
gpoints = 4
for i in range(0, nrNodes-1):
    test = zeros(4)
    points, weights = Gaussian(gpoints)
    for ii in range(0, gpoints):
        xi = lx * (1 + points[ii]) / 2
        nv, dnv, ddnv = Nv(xi, lx)
        test = test + (nv*load(lx*i+xi)) * lx * 0.5 * weights[ii]

    li = array([0, test[0], test[1]])
    lj = array([0, test[2], test[3]])

    p5.addLoad(1, i+1, li)
    p5.addLoad(1, i+2, lj)

    print(test)



p5.updateElements()


p5.partitionK()




target = array([0, array([0, -80, 0])])

start_time = time.time()
upath, load = p5.solveArc(load=10, s_lim=1000, nodeDis=target, la=.8, alfa=0.0)
print("--- %s seconds ---" % (time.time() - start_time))
#print(upath)
#print(load)
p5.plotStructure(deformed=True, mark=True)

start_time = time.time()

p5.partitionK()
#print(p5.kff)
#print(p5.rtot)
print("--- %s seconds ---" % (time.time() - start_time))



nrEle = nrNodes -1
mid = int(nrEle/2 * 3 + 1)
quart_1 = int(nrEle/4 * 3 + 1)
quart_3 = int(nrEle*3/4 * 3 + 1)

print(f'mid {mid}, quart {quart_1}, 3 quart {quart_3}')

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter(upath[:, mid], load, color='b', label='$at 1/2 u_v$ ')
ax1.scatter(upath[:, quart_1], load, color='r', label='$at 1/4 u_v$')
ax1.scatter(upath[:, quart_3], load, color='g', label='$at 3/4 u_v$')
plt.title(f'Homework 5-3, {nrEle} elements')
ax1.set_xlabel('Displacement [cm]')
ax1.set_ylabel('Load factor')
ax1.legend()

plt.show()
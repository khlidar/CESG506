from Structure import *



p6 = struct()

nrNodes = 33  # Number of nodes
length = 500  # Length between
h0 = length/100  # height at mid span
x = linspace(0, length, nrNodes)
y = (4*h0*(x/length))*(1-(x/length))  # curvature of beam
y = zeros_like(x)
for i in range(0, nrNodes):
    if i == 0:
        p6.addNode(i + 1, array([0, 0, 1]), array([x[i], y[i]]))
    elif i == nrNodes-1:
        p6.addNode(i + 1, array([1, 0, 1]), array([x[i], y[i]]))
    else:
        p6.addNode(i + 1, array([1, 1, 1]), array([x[i], y[i]]))

para = {'E': 1.0,
        'A': 10000.0,
        'I': 100000.0,
        'nu': 0.0}

for i in range(0, nrNodes - 1):
    p6.addEle(i + 1, i + 1, i + 2, CurvedBeam_V2(params=para))

# Assign loads
p6.addLoad(1, nrNodes, array([-1., 0., 0.]))

p6.updateElements()


p6.partitionK()


target = array([nrNodes, array([-100, 0, 0])])

start_time = time.time()
upath, load, s_hist, R_hist, det_hist = p6.solveArc(load=6, s_lim=500,nodeDis=target, la=0.2, alfa=0.0, data=True)
print("--- %s seconds ---" % (time.time() - start_time))

p6.plotStructure(deformed=True, mark=True)

start_time = time.time()

p6.partitionK()
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
#ax1.scatter(load[1::], det_hist[1::], color='b', label='$det[K]$ ')
ax1.plot(load[1::], det_hist[1::], color='b', label='$min(eig([Kff]))$ ')
ax1.plot([0, load[len(load)-1]], [0, 0], color='r')
ax1.plot([3.9478417604, 3.9478417604], [min(det_hist), max(det_hist)], label='True value straight beam')
#ax1.semilogy(load[1::], det_hist[1::], color='b', label='$det[K]$ ')
plt.title(f'Homework 6-4, {nrEle} elements straight beam')
ax1.set_xlabel('Load factor')
ax1.set_ylabel('Value')
#ax1.set(ylim=(0, 1000))
ax1.legend()

plt.show()
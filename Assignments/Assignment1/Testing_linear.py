


from numpy import array, linalg, outer, identity

k = array([[890, -29], [-29, 11]])
P = array([0., -0.1])

u = P/k

print(u)

P = array([[0.], [-0.1]])

u = P/k

print(u)

l = array([5.5, 0.5])

print(outer(l, l))

print(identity(3))

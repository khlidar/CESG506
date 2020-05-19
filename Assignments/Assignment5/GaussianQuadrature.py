
from scipy import integrate
from numpy import sqrt, array
import sys


def Gaussian(gpoints=2):
    if gpoints == 1:
        points = [0]
        weights = [2]
    elif gpoints == 2:
        points = [-0.5773502692, 0.5773502692]
        weights = [1, 1]
    elif gpoints == 3:
        points = [-0.7745966692, 0, 0.7745966692]
        weights = [0.5555555556, 0.8888888889, 0.5555555556]
    elif gpoints == 4:
        points = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116]
        weights = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]
    elif gpoints == 5:
        points = [-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459]
        weights = [0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851]
    else:
        print(f'Gaussian is only defined for 1-5, 5 points used instead of {gpoints}')
        points, weights = Gaussian(5)

    return points, weights


def k_test(x):
    j11 = 1/x
    j21 = -1/x

    j = array([[j11, j21],[j21, j11]])

    return j

def main():
    I = integrate.quad_vec(k_test, 1, 10)



    print(I)






if __name__ == "__main__":
    main()
    sys.exit(0)

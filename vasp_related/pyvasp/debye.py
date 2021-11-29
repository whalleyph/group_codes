import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

kb = 8.617343e-5        # eV/K

def integrand(x):
    exp = np.exp
    expx = exp(x)
    tmp = x**4*expx/(expx-1)**2
    return tmp

def integral(func, b):
    return quad(func, 0, b)

def cv(x):
    int = integral(integrand, 1/x)[0]
    int = 9*kb*x**3*int
    return int

def cvby3k(x):
    int = integral(integrand, 1/x)[0]
    int = 3*x**3*int
    return int

def cv_by_3k(temperature, td):
    cv = []
    for t in temperature:
        t_by_td = t/float(td)
        int = integral(integrand, 1/t_by_td)[0]
        int = 3*(t_by_td)**3*int
        cv.append(int)
    cv = np.array(cv)

    return cv

if __name__ == '__main__':

    tmin = 1
    tmax = 800
    temperature = np.arange(tmin, tmax, 1.0)

    for td in [474, 604]:
        cvby3k = cv_by_3k(temperature, td)
        plt.plot(temperature, cvby3k*3*kb)
        # print 3*kb*temperature[-1]*cvby3k[-1]
        print np.trapz(cvby3k[0:750]*3*kb, temperature[0:750])
    plt.legend(['bcc-Fe, 474 K', 'cementite, 604 K'], loc='lower right')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Cv [eV/K]')
    plt.savefig('Cv_bccFe_cementite.eps')
    # plt.show()

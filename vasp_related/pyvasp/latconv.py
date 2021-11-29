#! /usr/bin/python

from numpy.linalg import norm
import numpy as np

acos = np.arccos
cos = np.cos
pi = np.pi

def d2r(d):
    return pi*d/180.

def latvec2latpar(A):
    # A = [ax ay az, bx by bz, cx cy cz]
    av = A[0,:]
    bv = A[1,:]
    cv = A[2,:]

    a = norm(av)
    b = norm(bv)
    c = norm(cv)

    alpha = acos(np.dot(bv,cv)/(b*c))*180./pi
    beta  = acos(np.dot(av,cv)/(a*c))*180./pi
    gamma = acos(np.dot(av,bv)/(a*b))*180./pi

    return (a, b, c, alpha, beta, gamma)

def latpar2latvec(X):
    # X = (a, b, c, alpha, beta, gamma) 

    a = X[0]
    b = X[1]
    c = X[2]
    alpha = X[3]
    beta  = X[4]
    gamma = X[5]

    alpha = d2r(alpha)
    beta  = d2r(beta)
    gamma = d2r(gamma)

    ax = a
    ay = 0
    az = 0

    bz = 0

    bx = b*cos(gamma)
    by = (b**2-bx**2)**0.5

    cx = c*cos(beta)
    cy = (b*c*cos(alpha)-bx*cx)/by
    cz = (c**2-cx**2-cy**2)**.5

    return np.array([[ax, ay, az], [bx, by, bz], [cx, cy, cz]])




if __name__ == '__main__':
    A = np.array([[0, 0, 6.38],[6.419, 0, 0], [3.2, 5.67, 0]])
    latpar = latvec2latpar(A)
    print latpar
    print latpar2latvec(latpar)

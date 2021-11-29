def debye3(x):
    y = None
    if x <= 1.0:
        y = -4.635519*x*x - 2.333748e-3*x - 1.000327
    elif x < 4.0:
        y =  6.482976e-3*x*x*x \
            -5.684513e-2*x*x \
            -1.473625e-3*x \
            +1.003855
    elif x <=6.0:
        y =  1.566803e-2*x*x \
            -2.754144e-1*x \
            +1.354056
    else:
        print 'Function not defined for:', x

    return y

if __name__ == '__main__':
    print debye3(0.9)

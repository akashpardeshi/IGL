import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def momentK(G, k):
    '''
    Computes the kth moment function of the prime gaps
    param: G: a list of the first N prime gaps
    param: k: the moment to compute
    return: a list of the values for the function
    A_k(n) = 1/n sum(i = 1, n) G[i]^k
    '''
    N = len(G)
    Ak = N * [1]
    for i in range(1, N):
        Ak[i] = (Ak[i-1] * (i-1) + G[i]**k) / i
    return Ak

def f(x, a, b):
    # the exponent should match the k
    # in plotMomentK
    return a * np.log(x)**1 + b

def plotMomentK(G, k):
    '''
    plots the function A_k with the regression
    sampling from the gaps data
    param: G: a list of the first N prime gaps
    param: k: the moment
    '''
    N = len(G)
    # sample 1 in every 1000 gaps to
    # steed up plotting times
    step = 1000
    x = [i for i in range(1, N + 1, step)]
    G = [G[i] for i in range(N) if i % step == 0]
    y = momentK(G, k)
    assert(len(x) == len(y))

    # a, b = curve_fit(f, x, y)[0]
    # print(a, b)
    # assert(False)
    # parameters from regression for k = 1:
    a, b = 1.0649164354465381, 0.7746522519415989 
    # parameters from regression for k = 2:
    # a, b = 2.130410641914967, 9.332692023782698
    # parameters from regression for k = 3:
    # a, b = 6.159680892723837, -210.14343426983513
    # generate regression values for y
    yFit = [f(x, a, b) for x in x]

    # plotting and labelling
    plt.title(rf'$A_{k}$ mean function on the prime gaps', fontsize=20)
    plt.xlabel(r'$n$')
    plt.ylabel(rf'$A_{k}(n)$')
    plt.plot(x, y, 'o', alpha=0.25, markersize=5)
    plt.plot(x, yFit, 'k', label=rf'$y = a\log(n) + b$' '\n' rf'$a = {a:.6}, b = {b:.6}$')
    leg = plt.legend(loc='lower right', fontsize='large')
    leg.get_frame().set_linewidth(0.0)
    plt.show()

def main():
    G = []
    with open('gaps.txt') as gaps:
        for p in gaps:
            G.append(int(p.strip()))
    plotMomentK(G, 1)

main()

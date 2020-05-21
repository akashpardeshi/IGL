import matplotlib.pyplot as plt
import numpy as np
# import math
from scipy.optimize import curve_fit
import scipy.stats as s
import scipy.special as sps

plt.style.use('seaborn')

def collatz(N, found_seqs):
    if N == 1: # base case
        return 0
    if N not in found_seqs:
        found_seqs[N] = 1 + collatz(N//2 if N % 2 == 0 else 3*N + 1, found_seqs) # recursive step
    return found_seqs[N]

def aveLen(N):
    """
    returns a list of the average Collatz sequence lengths from 1 to n when n varies from 1 to N
    """
    found_seqs = {} # dictionary of number : Collatz length of number
    for n in range(1, N + 1): # populate the dictionary
        collatz(n, found_seqs)
    aveLengths = [0]
    for n in sorted(found_seqs): # calculate average lengths
        if n <= N:
            # aveLengths.append( (aveLengths[-1] * len(aveLengths) + found_seqs[n]) / (len(aveLengths)+1) )
            aveLengths.append( (aveLengths[-1] * (n-1) + found_seqs[n]) / n )
    return aveLengths

def maxNum(N):
    """
    returns the largest number in the collatz sequence of N
    """
    found_seqs = {} # dictionary of number : Collatz length of number
    collatz(N, found_seqs)
    return max(found_seqs)

def f(x, a, c):
    return a*np.log(x) + c

def g(x, a, b, c):
    return a*np.log(x)**c + b

def plotAveLen(N):
    """
    plots the average collatz length from 1...n when n varies from 1 to N
    """
    av_plt = plt.subplot()
    plt.title('Average Length of Collatz Sequences')
    plt.xlabel('N')
    plt.ylabel('average sequence length')

    n_vals = [n for n in range(1, N + 1)]
    aveLengths = aveLen(N)
    av_plt.plot(n_vals, aveLengths, 'o', alpha=0.5, label='scatter plot of average lengths')

    params = curve_fit(f, n_vals, aveLengths)
    a, c = params[0]
    print(a, c)
    y = [f(x, a, c) for x in n_vals]
    av_plt.plot(n_vals, y, '-k', label=f'curve fit following Alog(x) + C for A = {a:.3} and C = {c:.3}')
    av_plt.legend(loc='lower right')

    plt.show()

    return a, c

def plotAveResid(N):
    n_vals = [n for n in range(1, N + 1)]
    aveLengths = aveLen(N)

    res_plt = plt.subplot()
    plt.title('Residual Plot for Average Length of Collatz Sequences')
    plt.xlabel('N')
    plt.ylabel('Residuals')

    a, c = curve_fit(f, n_vals, aveLengths)[0]
    residuals = [aveLengths[n[0]] - f(n[1], a, c) for n in enumerate(n_vals)]
    res_plt.plot(n_vals, residuals, 'or', alpha=0.5, label='residuals for curve model')
    res_plt.legend(loc='upper right')

    plt.show()

def linReg(x, y, N):
    """
    returns a 2 lists of coefficients A, C for all n in 2, ..., N
    """
    t1 = time.time()
    a_coef = []
    c_coef = []
    sum_x = sum_y = sum_xy = sum_x2 = 0
    for n in range(1, N):
        sum_x += math.log(1+n)
        sum_y += y[n]
        sum_xy += math.log(1+n) * y[n]
        sum_x2 += math.log(1+n)**2
        a = ((1+n)*sum_xy - sum_x*sum_y) / ((1+n)*sum_x2 - sum_x**2)
        a_coef.append(a)
        c_coef.append( (sum_y - a*sum_x) / (n+1) )
    t2 = time.time()
    print(a_coef[-1], c_coef[-1])
    print(t2 - t1)
    return a_coef, c_coef

def plotParams(N):
    a_plt = plt.subplot(1, 2, 1)
    plt.title('Convergence of Parameter A')
    plt.xlabel('N')
    plt.ylabel('Value of A')

    n_vals = [n for n in range(1, N+1)]
    aveLengths = aveLen(N)
    a_vals, c_vals = linReg(n_vals, aveLengths, N)

    a_plt.plot(n_vals[1:], a_vals, 'm')
    a_plt.plot(n_vals, N*[10.348], 'k', linestyle='--')

    c_plt = plt.subplot(1, 2, 2)
    plt.title('Convergence of Parameter C')
    plt.xlabel('N')
    plt.ylabel('Value of C')

    c_plt.plot(n_vals[1:], c_vals, 'm')
    c_plt.plot(n_vals, N*[-11.552], 'k', linestyle='--')

    plt.show()

def plotMaxNum(N):
    """
    plots the largest number reached in the collatz sequence of N
    """
    plt.title('Largest Number Reached in Collatz Sequences')
    plt.xlabel('N')
    plt.ylabel('largest number reached')
    plt.yscale('log')
    n_vals = [n for n in range(1, N + 1)]
    maxNums = [maxNum(n) for n in n_vals]
    plt.plot(n_vals, maxNums, 'ro', alpha=0.1)
    plt.show()

def plotSeq(N):
    """
    plots the collatz length of N
    """
    found_seqs = {} # dictionary of number : Collatz length of number
    plt.title('Collatz Sequence Length')
    plt.xlabel('N')
    plt.ylabel('length')
    n_vals = [n for n in range(2, N + 1)]
    nums = [collatz(n, found_seqs) for n in n_vals]
    plt.plot(n_vals, nums, 'ro', alpha=0.25)
    plt.show()

def var_len(N):
    found_seqs = {}
    for n in range(1, N+1):
        collatz(n, found_seqs)
    lens = [0]
    for n in sorted(found_seqs):
        if n <= N:
            lens.append(found_seqs[n])
    var = []
    mean = aveLen(N)
    sum_x2 = 0
    for n in range(1, N+1):
        sum_x2 += lens[n-1]**2
        v = sum_x2/n - mean[n-1]**2
        var.append(v)
    return var

def plot_var(N):
    var = var_len(N)
    x = [x for x in range(1, N+1)]
    a, c = curve_fit(f, x, var)[0]
    y = [f(x, a, c) for x in x]
    plt.title('Variance of Collatz Lengths')
    plt.xlabel('x')
    plt.ylabel('variance')
    plt.plot(x, var, 'ko', alpha=0.5)
    plt.plot(x, y)
    plt.show()

    return a, c

def hist(N):
    found_seqs = {}
    for x in range(1, N+1):
        collatz(x, found_seqs) # populate the dictionary
    vals = [0] + [found_seqs[n] for n in found_seqs if n <= N]

    plt.title('Histogram of Collatz Sequence Lengths')
    plt.xlabel('Sequence Length')
    plt.ylabel('Frequency')

    x = np.linspace(1, max(vals), N)
    count, bins, ignored = plt.hist(vals, bins=100, color='purple', density=True)
    
    k, g, theta = s.gamma.fit(vals)
    y = s.gamma.pdf(bins, k, g, theta)
    plt.plot(bins, y, color='black', label=rf'Gamma distribution with k: {k:.5}, $\theta$: {theta:.5}')
    plt.legend(loc='upper right')

    plt.show()

def moments(N):
    # found_seqs = {} # dictionary of number : Collatz length of number
    # lens = [n for n in range(2, N + 1)]
    # nums = [collatz(n, found_seqs) for n in lens]
    lens = aveLen(N)
    m1 = [0.5]
    m2 = [0.5]
    m3 = [0.5]
    m4 = [0.5]
    m5 = [0.5]
    for n in range(1, len(lens)):
        m1.append((m1[-1]*n + lens[n])/(n+1))
        m2.append((m2[-1]*n + lens[n]**2)/(n+1))
        m3.append((m3[-1]*n + lens[n]**3)/(n+1))
        m4.append((m4[-1]*n + lens[n]**4)/(n+1))
        m5.append((m5[-1]*n + lens[n]**5)/(n+1))
    # print(m1, m2, m3, m4, m5)
    return m1, m2, m3, m4, m5

def plot_moments(N):
    m1, m2, m3, m4, m5 = moments(N)
    x = [x for x in range(1, len(m1) + 1)]
    # plt.plot(x, m1)
    # plt.plot(x, m2)
    # plt.plot(x, m3)
    # plt.plot(x, m4)
    # plt.plot(x, m5)
    # plt.show()
    a1, b1, c1 = curve_fit(g, x, m1)[0]
    a2, b2, c2 = curve_fit(g, x, m2)[0]
    a3, b3, c3 = curve_fit(g, x, m3)[0]
    a4, b4, c4 = curve_fit(g, x, m4)[0]
    a5, b5, c5 = curve_fit(g, x, m5)[0]
    print(f'moment 1: {a1, b1, c1}')
    print(f'moment 2: {a2, b2, c2}')
    print(f'moment 3: {a3, b3, c3}')
    print(f'moment 4: {a4, b4, c4}')
    print(f'moment 5: {a5, b5, c5}')

def main():
    N = int(input("Enter the upper bound for N: "))
    # plotAveResid(N)
    # linReg(n_vals, aveLengths, N)
    # plotParams(N)
    # plotAveLen(N)
    # plotMaxNum(N)
    # plotSeq(N)
    # plot_var(N)
    # hist(N)
    plot_moments(N)

main()

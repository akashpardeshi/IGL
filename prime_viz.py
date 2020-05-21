import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as s


# plt.style.use('seaborn-poster')
plt.style.use('seaborn')
# plt.style.use('ggplot')

def sieve(N):
    """
    returns a list of where the primes are under N (1 indicates primes, 0 indicates composite)
    """
    nums = (N+1) * [1]
    for n in range(2, int(N**0.5) + 1):
        if nums[n] == 0:
            continue
        k = n
        while k + n <= N: # mark off the multiples of n
            k += n
            nums[k] = 0
    nums[0] = nums[1] = 0
    return nums

def sieve2(N):
    """
    returns a list of where the primes are under N (1 indicates primes, 0 indicates composite)
    """
    nums = [x for x in range(N+1)]
    for n in range(2, int(N**0.5) + 1):
        if nums[n] == 0:
            continue
        k = n
        while k + n <= N: # mark off the multiples of n
            k += n
            nums[k] = 0
    nums[0] = nums[1] = 0
    return [x for x in nums[2:] if x != 0]
 
def pi(N):
    p = sieve(N)
    y = [0]
    for x in p[1:]:
        y.append(y[-1]) if x == 0 else y.append(y[-1] + 1)
    return y

def pnt(N):
    x = [x for x in range(N + 1)]
    y = pi(N)
    a = [x / np.log(x) for x in x[2:]]

    pnt1 = plt.subplot(1, 2, 1)
    plt.title('Prime Number Theorem', fontsize=20)
    plt.xlabel(r'$x$', fontsize=15)
    plt.ylabel(r'number of primes less than or equal to $x$', fontsize=15)

    plt.step(x, y, 'c', linewidth=2, where='post', label=r'$\pi(x)$')
    plt.plot(x[2:], a, 'r', linewidth=2, label=r'$\frac{x}{\logx}$')
    leg = plt.legend(loc='lower right', fontsize='x-large')
    leg.get_frame().set_linewidth(0.0)

    pnt2 = plt.subplot(1, 2, 2)
    plt.ylim(0.95, 1.3)
    plt.title(r'Asymptotic relationship' '\n' r'between $\pi(x)$ and $\frac{x}{\logx}$', fontsize=20)
    plt.xlabel(r'$x$', fontsize=15)
    plt.ylabel(r'ratio of $\pi(x)$ to $\frac{x}{\logx}$', fontsize=15)

    plt.plot(x[2:], [y/a for y, a in zip(y[2:], a)], 'y', linewidth=2)

    plt.show()

def gap(N):
    primes = [p[0] for p in enumerate(sieve(N)) if p[1] != 0]
    # primes = [p[0] for p in enumerate(sieve2(N)) if p[1] != 0]
    gaps = [primes[i+1] - primes[i] for i in range(len(primes) - 1)]
    return gaps

def plot_gap(N):
    gaps = gap(N)
    x = [x for x in range(1, len(gaps) + 1)]

    plt.title(r'Prime Gaps', fontsize=20)
    plt.xlabel(r'$n$', fontsize=15)
    plt.ylabel(r'$\gamma(n)$', fontsize=15)

    plt.plot(x, gaps, 'g.', alpha=0.35, markersize=2)
    plt.show()

def ave_gap(N):
    """
    returns a list of the averge prime gap to the Nth prime
    """
    gaps = gap(N)
    ave_gaps = [1]
    for n in range(1, len(gaps)):
        ave_gaps.append( (ave_gaps[-1]*n + gaps[n]) / (n + 1) )
    return ave_gaps

    # t1 = time.time()
    # primes = [p[0] for p in enumerate(sieve(N)) if p[1] != 0]
    # ave_gaps = [(primes[i] - primes[0])/i for i in range(1, len(primes))]
    # t2 = time.time()
    # print(t2 - t1)
    # return ave_gaps

def g(x, a, b, c):
    return a*(np.log(1+x)**c) + b

def f(x, a, b):
    return a*np.log(1+x) + b

def plot_ave_gap(N):
    ave_gaps = ave_gap(N)
    x = [x for x in range(len(ave_gaps))]
    # a, b, c = curve_fit(f, x, ave_gaps, p0=[1, 1, 1.4])[0]
    # print(a, b, c)
    # y = [f(x, a, b, c) for x in x]

    # plt.subplot(1, 2, 1)
    plt.title('Cumulative Mean of Prime Gaps', fontsize=20)
    plt.xlabel(r'$n$', fontsize=15)
    plt.ylabel(r'$\alpha(n)$', fontsize=15)
    plt.plot(x, ave_gaps, 'ko', alpha=0.5, markersize=5)
    # plt.plot(x, y, 'b')

    # plt.subplot(1, 2, 2)
    # ratios = [y/g for y, g in zip(y, ave_gaps)]
    # plt.plot(x, ratios)
    # plt.plot(x, len(x)*[1], 'k--')

    plt.show()
    # return a, b, c

def plot_mean_var(N):
    ave_gaps = ave_gap(N)
    # var = var_gap(N)
    x = [x for x in range(len(ave_gaps))]

    fig1 = plt.subplot(1, 2, 1)
    a, b, c = curve_fit(g, x, ave_gaps)[0]
    print(a, b, c)
    y = [g(x, a, b, c) for x in x]
    plt.title('Cumulative Mean of Prime Gaps', fontsize=20)
    plt.xlabel(r'$n$', fontsize=15)
    plt.ylabel('Mean', fontsize=15)
    plt.plot(x, ave_gaps, '.', markersize=10, color="C1", label=r'$\alpha(n)$') # plot data
    plt.plot(x, y, 'k', linewidth=2, label=rf'$y = a\log(n)^c + b$' '\n' rf'$a = {a:.4}, b = {b:.4}, c = {c:.4}$') # plot regression
    leg = plt.legend(loc='lower right', fontsize='x-large')
    leg.get_frame().set_linewidth(0.0)

    plt.subplot(1, 2, 2)
    plt.ylim(-0.2, 0.2)
    plt.title(r'Difference Between $\alpha(n)$' '\n' 'and Regression Model', fontsize=20)
    plt.xlabel(r'$n$', fontsize=15)
    plt.ylabel('Residuals', fontsize=15)
    difs = [a-y for y, a in zip(y, ave_gaps)]
    plt.plot(x, difs, '.', alpha=0.1, markersize=5, label=r'$\alpha(n) - (a\log(n)^c + b)$')
    plt.plot(x, len(x)*[0], 'k--', markersize=5)
    leg = plt.legend(loc='lower right', fontsize='large')
    leg.get_frame().set_linewidth(0.0)

    plt.show(fig1)
#############
    # fig2 = plt.subplot(1, 2, 1)
    # plt.subplot(1, 2, 2)
    # # a, b, c = curve_fit(g, x, var)[0]
    # # print(a, b, c)
    # # y = [g(x, a, b, c) for x in x]
    # plt.title('Cumulative Variance of Prime Gaps', fontsize=20)
    # plt.xlabel(r'$n$', fontsize=15)
    # plt.ylabel('Variance', fontsize=15)
    # plt.plot(x, var, '.', markersize=10, label=r'$\beta(n)$') # plot data
    # # plt.plot(x, y, 'k', markersize=5, label=rf'$y = a\log(n)^c + b$' '\n' rf'$a = {a:.6}, b = {b:.6}, c = {c:.6}$') # plot regression
    # leg = plt.legend(loc='lower right', fontsize='x-large')
    # leg.get_frame().set_linewidth(0.0)

    # # plt.subplot(1, 2, 2)
    # # ratios = [v-y for y, v in zip(y, var)]
    # # plt.title(r'Ratio of $\beta(n)$ to Regression Model', fontsize=20)
    # # plt.xlabel(r'$n$', fontsize=15)
    # # plt.ylabel(r'$\frac{\beta(n)}{a\log(n)^c + b}$', fontsize=15)
    # # plt.plot(x, ratios, 'ro', alpha=0.5, markersize=5)
    # # plt.plot(x, len(x)*[1], 'k--')

    # # plt.show(fig2)
    # plt.show()

def var_gap(N):
    mean = ave_gap(N)
    gaps = gap(N)
    var = []
    sum_x2 = 0
    for n in range(1, len(mean) + 1):
        sum_x2 += gaps[n-1]**2
        v = sum_x2/n - mean[n-1]**2
        var.append(v)
    return var

def plot_var(N):
    var = var_gap(N)
    x = [x for x in range(1, len(var) + 1)]
    a, b, c = curve_fit(g, x, var)[0]
    print(a, b, c)
    y = [g(x, a, b, c) for x in x]

    plt.subplot(1, 2, 1)
    # plt.plot(x[50000:], var[50000:], 'ko')
    plt.plot(x, var, 'ko', alpha=0.5, markersize=5)
    plt.plot(x, y)
    plt.title('Cumulative Variance of Prime Gaps', fontsize=20)
    plt.xlabel(r'$n$', fontsize=15)
    plt.ylabel(r'$\beta(n)$', fontsize=15)

    plt.subplot(1, 2, 2)
    ratios = [v/y for y, v in zip(y, var)]
    plt.plot(x, ratios, 'o')
    plt.plot(x, len(x)*[1], 'k--')

    plt.show()
    # return a, b, c

def hist(N):
    vals = gap(N)[1:]
    vals = [v/2 for v in vals]
    plt.title('Histogram of Prime Gaps')
    plt.ylabel('Frequency')

    mu = np.mean(vals)
    sig = np.std(vals)
    n, bins, patches = plt.hist(vals, bins=100, density=True)

    a, b, c = s.lognorm.fit(vals)
    y = s.lognorm.pdf(bins, a, b, c)
    plt.plot(bins, y)

    plt.show()

def moments(N):
    gaps = gap(N)
    m1 = [1]
    m2 = [1]
    m3 = [1]
    m4 = [1]
    m5 = [1]
    for n in range(1, len(gaps)):
        m1.append((m1[-1]*n + gaps[n])/(n+1))
        m2.append((m2[-1]*n + gaps[n]**2)/(n+1))
        m3.append((m3[-1]*n + gaps[n]**3)/(n+1))
        m4.append((m4[-1]*n + gaps[n]**4)/(n+1))
        m5.append((m5[-1]*n + gaps[n]**5)/(n+1))
    return m1, m2, m3, m4, m5

def hi(N):
    m1, m2, m3, m4, m5 = moments(N)
    x = [x for x in range(1, len(m1) + 1)]
    a1, b1, c1 = curve_fit(g, x, m1)[0]
    a2, b2, c2 = curve_fit(g, x, m2)[0]
    a3, b3, c3 = curve_fit(g, x, m3)[0]
    a4, b4, c4 = curve_fit(g, x, m4)[0]
    a5, b5, c5 = curve_fit(g, x, m5)[0]
    print(f'moment 1: {c1}')
    print(f'moment 2: {c2}')
    print(f'moment 3: {c3}')
    print(f'moment 4: {c4}')
    print(f'moment 5: {c5}')
    plt.show()

def error(N):
    a, b, c = 1.2785858207407403, -0.12907491094053686, 0.953059290922383 
    primes = sieve2(N)
    x = [x for x in range(1, len(primes) + 1)]
    # y1 = np.abs([primes[n] - (n+1)*np.log(n+1) for n in range(1, len(x))])
    y2 = np.abs([primes[n] - (primes[n-1] + g(n-1, a, b, c)) for n in range(1, len(x))])

    # fig1 = plt.title('Errors in Predictions of Primes with Analytical and Probabilistic Approaches', fontsize=20)
    # label = plt.xlabel(r'$n$', fontsize=15)
    # plt.ylabel(r'$\vert p_{n+1} - \hat{p}_{n+1} \vert$', fontsize=15)

    # plt.plot(x[:-1], y1, '.', color='orchid', alpha=0.5, markersize=5, label=r'Errors when $\hat{p}_{n+1} = (n+1)\log(n+1)$')
    # plt.plot(x[:-1], y2, '.', color='lightseagreen', alpha=0.5, markersize=5, label=r'Errors when $\hat{p}_{n+1} = p_{n} + \alpha(n)$')
    # leg = plt.legend(loc='upper left', fontsize='large')
    # leg.get_frame().set_linewidth(0.0)
    # plt.show(fig1)

    # fig2 = plt.subplot(1, 2, 1)
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    # plt.title('Histogram of Errors in Analytic \n Predictions of Primes', fontsize=20)
    # plt.xlabel(r'Errors: $\vert p_{n+1} - \hat{p}_{n+1} \vert$', fontsize=15)
    # plt.ylabel('Frequency', fontsize=15)
    # n, bins, patches = plt.hist(y1, bins=50, color='orchid', label=r'Errors when $\hat{p}_{n+1} = (n+1)\log(n+1)$')
    # plt.xticks(np.arange(0, max(bins)+1, 200000))
    # leg = plt.legend(fontsize='large')
    # leg.get_frame().set_linewidth(0.0)
    # # plt.show(fig2)

    # plt.subplot(1, 2, 2)
    plt.title('Histogram of Errors in Probabilistic \n Predictions of Primes', fontsize=20)
    plt.xlabel(r'Errors: $\vert p_{n+1} - \hat{p}_{n+1} \vert$', fontsize=15)
    plt.ylabel('Frequency', fontsize=15)
    n, bins, patches = plt.hist(y2, bins=50, color='lightseagreen', label=r'Errors when $\hat{p}_{n+1} = p_{n} + \alpha(n)$')
    plt.xticks(np.arange(0, max(bins)+1, 20))
    leg = plt.legend(fontsize='large')
    leg.get_frame().set_linewidth(0.0)
    plt.show()

def asy(N):
    r1 = ave_gap(N)[-1]
    primes = sieve2(N)
    r2 = sum([np.log(p) for p in primes]) / len(primes)
    return r1, r2

def plot(N):
    x = np.arange(2, N)
    y = [] 
    for n in range(2, N):
        r1, r2 = asy(n)
        y.append(r1/r2)
    plt.plot(x, y)
    plt.show()

def main():
    N = int(input('Enter a value for N: '))
    # pnt(N)
    # gap(N)
    # plot_gap(N)
    # ave_gap(N)
    # plot_ave_gap(N)
    # print(ave_gap(N))
    # plot_var(N)
    # plot_mean_var(N)
    # hist(N)
    # error(N)
    # hi(N)
    # asy(N)
    plot(N)


main()

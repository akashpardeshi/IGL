import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.style.use('ggplot')

######################################################
##### Functions for generating prime number data #####
######################################################

def generatePrimes(N, bit_vector=False):
    """
    Returns a list of all primes <= N using sieve of Eratosthenes.
    bit_vector: False results in a list of prime numbers
                True results in a list where list[n] == 1 if n is prime else 0
    """
    is_prime = [0, 0, 1] + [1 if n % 2 else 0 for n in range(3, N + 1)]
    for p in range(3, int(N**0.5) + 1):
        if is_prime[p]:
            for m in range(p**2, N + 1, p):
                is_prime[m] = 0
    return is_prime if bit_vector else [n for n in range(2, N + 1) if is_prime[n]]

def generateGaps(N):
    """
    Returns a list of the prime gaps for primes <= N.
    """
    primes = generatePrimes(N)
    return [primes[n+1] - primes[n] for n in range(len(primes) - 1)]

def pi(N):
    """
    Returns a list where the nth index is pi(n) := number of primes <= n for n = 0..N
    """
    primes = generatePrimes(N, bit_vector=True)
    pi = [0]
    for n in primes[1:]:
        pi.append(pi[-1] + n)
    return pi

def kthMomentPrimeGaps(N, k):
    """
    Returns a list where list[i] = kth moment of [prime gaps for primes <= i] for i = 2..N
    Define the kth moment of a sequence [a_1..a_n] := 1/n * sum(a_i**k for i = 1..n).
    """
    gaps = generateGaps(N)
    gaps_len = len(gaps)
    moments = [0 for _ in range(gaps_len)]
    for i in range(2, gaps_len):
        moments[i] = moments[i-1] + gaps[i]**k
    return [moments[i] / i for i in range(2, gaps_len)]

####################################################
##### Functions for plotting prime number data #####
####################################################

def plotPrimeNumberTheorem(N, save_fig=False):
    """
    Plots the pi(n) with its asymptotic relation to n / log(n).
    """
    pnt_data = np.array(pi(N))
    n_data = np.array(range(2, len(pnt_data)))
    log_data = n_data / np.log(n_data)
    pnt_data = pnt_data[2:]
    ratio_data = pnt_data / log_data

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
    fig.suptitle('Prime Number Theorem')

    ax1.set_title('Prime Counting Function')
    ax1.set(xlabel='$n$')
    ax1.step(n_data, pnt_data, 'c', where='post', label='$\pi(n)$')
    ax1.plot(n_data, log_data, 'r', label=r'$\frac{n}{\log n}$')
    ax1.legend(loc='lower right')

    ax2.set_title(r'Asymptotic Relationship of $\pi(n)$ and $\frac{n}{\log n}$')
    ax2.set(xlabel='$n$')
    ax2.plot(n_data, ratio_data, 'y', label=r'$\frac{\pi(n)}{\frac{n}{\log n}}$')
    ax2.legend(loc='lower right')

    if save_fig:
        plt.savefig('figures/prime_number_theorem.png')
    plt.show()

def plotPrimeGaps(N, save_fig=False):
    """
    Plots the prime gaps for primes <= N.
    """
    gap_data = generateGaps(N)
    
    plt.title('Prime Gaps')
    plt.xlabel('$n$')
    plt.ylabel('$n$th prime gap')
    plt.plot(range(1, len(gap_data) + 1), gap_data, 'g.', alpha=0.35, markersize=2)
    if save_fig:
        plt.savefig('figures/prime_gaps.png')
    plt.show()

##########################################################
##### Functions for plotting prime number statistics #####
##########################################################

def plotPrimeGapsMoment(N, k, save_fig=False):
    """
    Plots the kth moment data for each n = 1..N
    """
    moment_data = kthMomentPrimeGaps(N, k)
    moment_data_len = len(moment_data)
    n_data = np.array(range(1, moment_data_len+1))

    f = lambda x, a, b, c: a * np.log(x+1)**c + b
    params = curve_fit(f, n_data, moment_data)
    a, b, c = params[0]
    model_data = [f(n, a, b, c) for n in n_data]
    ratio_data = [d2 / d1 for d1, d2 in zip(model_data, moment_data)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))

    ax1.set_title(f'Moment-{k} for Primes')
    ax1.set(xlabel='$n$')
    ax1.plot(n_data, moment_data, '.', markersize=7, color='C1', label='moment data')
    ax1.plot(n_data, model_data, '-k', label=f'curve fit following $A\log(n+1)^C + B$\n(A = {a:.5}, B = {b:.5}, C = {c:.5})')
    ax1.legend(loc='lower right')

    ax2.set_title(f'Ratio of Model to Moment-{k} for Primes')
    ax2.set(xlabel='$n$', ylabel='ratio')
    ax2.plot(n_data, ratio_data, '.', markersize=7, color='C0')
    ax2.plot(n_data, [k for _ in range(moment_data_len)], '-k', label='1')
    ax2.legend(loc='upper right')

    if save_fig:
        plt.savefig(f'figures/prime_gaps_moment_{k}.png')
    plt.show()

def plotPrimeGapsHistogram(N, save_fig=False):
    """
    Plots histogram of prime gaps for primes <= N.
    """
    gaps = generateGaps(N)

    plt.title('Distribution of Prime Gaps')
    plt.xlabel('gap size')
    plt.ylabel('frequency')
    plt.hist(gaps, bins='fd', color='lightseagreen')
    if save_fig:
        plt.savefig('figures/prime_gaps_histogram.png')
    plt.show()

def main():
    N, k = 10**7, 2
    # plotPrimeNumberTheorem(N)
    # plotPrimeGaps(N)
    # plotPrimeGapsMoment(N, k)
    # plotPrimeGapsHistogram(N)

main()

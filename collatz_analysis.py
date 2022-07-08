import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.style.use('ggplot')

##########################################################
##### Functions for generating Collatz sequence data #####
##########################################################

def collatzSeqLen(N, cache):
    """
    Returns the length of the Collatz sequence starting at N.
    cache: dictionary for caching known collatz sequence lengths.
    """
    if N == 1:
        cache[N] = 1
    if N not in cache:
        cache[N] = 2 + collatzSeqLen((3*N + 1) // 2, cache) if N % 2 else 1 + collatzSeqLen(N // 2, cache)
    return cache[N]

def collatzSeq(N, cache):
    """
    Returns the Collatz sequence starting at N as a list.
    cache: dictionary for caching known Collatz sequences
    """
    if N == 1:
        cache[N] = [1]
    if N not in cache:
        cache[N] = [N] + [3*N + 1] + collatzSeq((3*N + 1) // 2, cache) if N % 2 else [N] + collatzSeq(N // 2, cache)
    return cache[N]

def kthMomentCollatz(N, k):
    """
    Returns a list where list[i] = kth moment of [collatzSeqLen(n) for n = 1..i] for i = 1..N
    Define the kth moment of a sequence [a_1..a_n] := 1/n * sum(a_i**k for i = 1..n).
    """
    collatz_len = {}
    for n in range(1, N+1):
        collatzSeqLen(n, collatz_len)

    moments = [0 for _ in range(N+1)]
    for i in range(1, N+1):
        moments[i] = (moments[i-1] * (i-1) + collatz_len[i]**k) / i
    return moments

########################################################
##### Functions for plotting Collatz sequence data #####
########################################################

def plotCollatzSeq(N, save_fig=False):
    """
    Plots the elements of the Collatz sequence starting at N.
    """
    collatz_seq = collatzSeq(N, {})

    plt.title(f'Collatz Sequence Starting at {N}')
    plt.xlabel('$n$')
    plt.ylabel('$n$th element of sequence')
    plt.yscale('log')
    plt.plot(range(1, len(collatz_seq)+1), collatz_seq, '.', markersize=7)
    if save_fig:
        plt.savefig('figures/collatz_seq.png')
    plt.show()

def plotCollatzSeqLength(N, save_fig=False):
    """
    Plots the length of the Collatz sequence starting at n for each n = 1..N
    """
    collatz_len = {}
    collatz_seq_lengths = [collatzSeqLen(n, collatz_len) for n in range(1, N+1)]

    plt.title('Collatz Sequence Lengths')
    plt.xlabel('$n$')
    plt.ylabel('length of sequence')
    plt.plot(range(1, N+1), collatz_seq_lengths, 'm.', markersize=5, alpha=0.1)
    if save_fig:
        plt.savefig('figures/collatz_seq_lengths.png')
    plt.show()

def plotMaxNumInCollatzSeq(N, save_fig=False):
    """
    Plots the largest number in Collatz sequence for each n = 1..N
    """
    collatz_seq = {}
    max_in_seq = [max(collatzSeq(n, collatz_seq)) for n in range(1, N+1)]

    plt.title('Largest Number Reached in Collatz Sequences')
    plt.xlabel('$n$')
    plt.ylabel('largest number reached')
    plt.yscale('log')
    plt.plot(range(1, N+1), max_in_seq, 'c.', markersize=5, alpha=0.1)
    if save_fig:
        plt.savefig('figures/max_num_in_collatz_seq.png')
    plt.show()

##############################################################
##### Functions for plotting Collatz sequence statistics #####
##############################################################

def plotCollatzSeqLengthMoment(N, k, save_fig=False):
    """
    Plots the kth moment data for each n = 1..N
    """
    moment_data = kthMomentCollatz(N, k)[1:]
    n_data = np.array(range(1, N+1))
    
    f = lambda x, a, b, c: a * np.log(x+1)**c + b
    params = curve_fit(f, n_data, moment_data)
    a, b, c = params[0]
    model_data = [f(n, a, b, c) for n in n_data]
    ratio_data = [d2 / d1 for d1, d2 in zip(model_data, moment_data)]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))

    ax1.set_title(f'Moment-{k} for Collatz Sequence Lengths')
    ax1.set(xlabel='$n$', ylabel='sequence length')
    ax1.plot(n_data, moment_data, '.', markersize=7, color='C1', label='moment data')
    ax1.plot(n_data, model_data, '-k', label=f'curve fit following $A\log(n+1)^C + B$\n(A = {a:.5}, B = {b:.5}, C = {c:.5})')
    ax1.legend(loc='lower right')

    ax2.set_title(f'Ratio of Model to Moment-{k} of Collatz Sequence Lengths')
    ax2.set(xlabel='$n$', ylabel='ratio')
    ax2.plot(n_data, ratio_data, '.', markersize=7, color='C0', label='moment data')
    ax2.plot(n_data, [1 for _ in range(N)], '-k', label='1')
    ax2.legend(loc='lower right')

    if save_fig:
        plt.savefig(f'figures/collatz_seq_length_moment{k}.png')
    plt.show()

def plotCollatzSeqLengthHistogram(N, save_fig=False):
    """
    Plots histogram of Collatz seqence lengths for sequences starting at 1..N
    """
    collatz_len = {}
    collatz_seq_lengths = [collatzSeqLen(n, collatz_len) for n in range(1, N+1)]

    plt.title('Distribution of Collatz Sequence Lengths')
    plt.xlabel('sequence length')
    plt.ylabel('frequency')
    plt.hist(collatz_seq_lengths, bins='fd', color='lightseagreen')
    if save_fig:
        plt.savefig(f'figures/collatz_seq_length_histogram.png')
    plt.show()

def main():
    N, k = 10**6, 1
    # plotCollatzSeq(N)
    # plotCollatzSeqLength(N)
    # plotMaxNumInCollatzSeq(N)
    # plotCollatzSeqLengthMoment(N, k)
    # plotCollatzSeqLengthHistogram(N)

main()

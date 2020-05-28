'''
generate a list of all primes <= N
and a list of the corresponding prime gaps
'''

import math

def generatePrimes(N):
    '''
    implementation of sieve of Eratosthenes
    return: a list of all primes <= N
    '''
    primes = [True for i in range(N+1)]
    for p in range(2, int(math.sqrt(N)) + 1):
        if p:
            for idx in range(p**2, N+1, p):
                primes[idx] = False

    return [p for p in range(2, len(primes)) if primes[p]]

def main():
    primes = generatePrimes(10**9)
    gaps = [primes[n+1] - primes[n] for n in range(len(primes) - 1)]
    with open('primes.txt', 'w') as data:
        for p in primes:
            data.write(f'{p}\n')
    with open('gaps.txt', 'w') as data:
        for g in gaps:
            data.write(f'{g}\n')

main()

def sieve(N):
    """
    returns a list of the primes less than or equal to N
    """
    nums = [n for n in range(N + 1)]
    for n in range(2, int(N**0.5) + 1):
        if nums[n] == 0:
            continue
        k = n
        while k + n <= N: # mark off the multiples of n
            k += n
            nums[k] = 0
    return [n for n in nums[2:] if nums[n] != 0]

def factor(N):
    """
    returns a list of the prime factors of N
    """
    factors = [] # list of prime factors of N
    primes = sieve(N) # a list of the candidates for prime factors
    k = 0
    while N > 1:
        for n in primes[k:]:
            if N % n == 0:
                factors.append(n)
                N /= n
                break
            k += 1 # we have found all the n's in the factorization of N
    return factors

def main():
    N = 2019
    # print(sieve(N))
    print(factor(N))

main()


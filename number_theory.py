"""
A helper library that will contain various
functions that are helpful in number theory.
Built for my convinience as I teach myself
number theory.
"""

def gcd(a, b):
    """
    Calculates the greatest common divisor between two numbers
    """
    if a == b and b == 0: return -1

    if a < b: return gcd(b, a)
    if b == 0: return a
    r = a % b
    return gcd(b, r)


def linear_gcd(a, b):
    """
    Computes gcd(a, b) as well as values x and y that satisfy ax + by = gcd(a, b)
    """
    if a < b: return linear_gcd(b, a)

    if b == 0:
        return a, 0, a

    x = 1
    g = a
    v = 0
    w = b
    while w != 0:
        q = g // w
        t = g % w
        s = x - q * v
        x = v
        g = w
        v = s
        w = t

    return g, (g - a * x) // b, x


def positive_linear_gcd(a, b):
    """
    Modified version of linear_gcd so that it always returns a solution with x > 0
    """
    g, x, y = linear_gcd(a, b)
    while x <= 0:
        x = x + b
        y = y - b
    return g, x, y


def factor(n):
    """
    Factors the number n and returns an list of lists containing the prime
    factors of n and the powers they are raised to
    """
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
		47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
		137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
		227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
		313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
		419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
		509, 521, 523, 541]
    factors = []
    for p in primes:
        if n % p != 0: continue
        n //= p
        factors.append([p, 1])
        while n % p == 0:
            n //= p
            factors[-1][1] += 1

    for d in range(primes[-1] + 2, n + 1, 2):
        if n % d != 0: continue
        n //= d
        factors.append([d, 1])

        while n % d == 0:
            n //= d
            factors[-1][1] += 1

    return factors


def print_factorization(n):
    """
    Factors the number n and prints out the prime factorization in a pleasing manner
    """
    factors = factor(n)
    for f in factors:
        print(f[0], end="")
        if f[1] != 1: print("^", f[1], sep="", end="")
        if f != factors[-1]: print(" * ", end="")


def solve_linear_congruence(a, c, m):
    """
    Solves linear congruences in the form of ax === c (mod m)
    """
    g, u, v = linear_gcd(a, m)
    if c % g != 0: return []

    solutions = []
    for i in range(g):
        x = u * (c // g)
        solutions.append((x + i * (m // g)) % m)

    return solutions


def solve_polynomial_congruence(coefficients, m):
    """
    Given the coefficients of the polynomial f(x), solves the congruence in the
    form f(x) === 0 (mod m)
    """
    solutions = []
    for i in range(m):
        value = 0
        for j in range(coefficients):
            value += (coefficients[j] * i ** (len(coefficients) - j - 1)) % m
        if value % m == 0: solutions.append(i)
    return solutions


def phi(n):
    """
    Computes Euler's phi function, counts the numbers less than or equal to n
    that are relatively prime to n
    """
    primes = factor(n)
    value = n;
    for p in primes:
        value *= 1 - 1 / p[0]
    return int(value)


def sigma(n):
    """
    Computes the sigma function, which is the sum of the divisors of n
    """
    factors = factor(n)
    value = 1
    for f in factors:
        value *= (f[0] ** (f[1] + 1) - 1) // (f[0] - 1)
    return value


def succ_square_power(a, k, m):
    """
    Computes a^k (mod m) using the method of successive squaring
    """
    b = 1
    while k >= 1:
        if k % 2 != 0: b = a * b % m
        a = a ** 2 % m
        k = k // 2
    return b


def root_mod(k, b, m):
    """
    Finds the roots of b modulo m, i.e solves the congruence in the form
    x^k === b (mod m)
    """
    g, u, v = positive_linear_gcd(k, phi(m))
    return succ_square_power(b, u, m)


def is_carmichael(n):
    """
    Determines if n is a Carmichael number
    """
    factors = factor(n)
    if len(factors) == 1: return False

    for f in factors:
        if f[1] != 1: return False
        if (n - 1) % (f[0] - 1) != 0: return False
    return True


def rabin_miller_test(n, threshold=10):
    """
    Uses the Rabin-Miller test to determine if the given number is composite
    Returns True if number is composite, and False if number is probably prime.
    The threshold value can be adjusted to check for more numbers to find witnesses
    """
    if n % 2 == 0: return True
    
    m = n - 1
    k = 0
    while m % 2 == 0:
        m //= 2
        k += 1
    q = (n - 1) // 2 ** k

    t = 0
    for a in range(2, n):
        if gcd(a, n) != 1: continue
        
        t += 1
        if t > threshold: break
        
        flag = False
        if succ_square_power(a, q, n) == 1: flag = True
        for i in range(0, k):
            if succ_square_power(a, 2 ** i * q, n) == n - 1: flag = True
        if flag: continue

        return True

    return False


print(rabin_miller_test(118901521, 10))

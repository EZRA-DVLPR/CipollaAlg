import argparse
from sympy import sqrt, expand, simplify

#returns a boolean on whether the given number is prime or not
def isPrime(n):
    #edge cases: n <= 1
    if n <= 1: return False

    #iterative approach
    for i in range(2, n):
        if n % i == 0:
            return False
        
    return True

#determines (2/p) and (3/p) given the Quadratic Natures of 2 and 3
def quadNat(num, p):
    if num != 2 and num != 3:
        raise ValueError ('Expected 2 or 3 but instead received ', num)
    
    elif num == 2:
        if 1 == p % 8 or 7 == p % 8:
            return 1
        else:
            #3 == p % 8 or 5 == p % 8
            return -1
    elif num == 3:
        if 1 == p % 12 or 11 == p % 12:
            return 1
        else: 
            #5 == p % 12 or 7 == p mod 12
            return -1

#determines the Legendre Symbol (t/p)
def legSymb(t, p):

    #if t is negative, then we separate the negative and return the Legendre symbol of the positive part
    if t < 0:
        #calculate using the Quadratic Nature of -1
        if 1 == p % 4:
            return legSymb(t * -1, p)
        else:
            return legSymb(t * -1, p) * (-1)

    #determine if t = 1, 2^k, or 3^k for some k and return the result of the Legendre Symbol
    if t == 1:
        return 1
    elif t % 2 == 0:
        return quadNat(2, p) * legSymb(int(t / 2), p)
    elif t % 3 == 0:
        return quadNat(3, p) * legSymb(int(t / 3), p)

    #t is not a power of 2 or 3, so perform the Law of Quadratic Reciprocity or Jacobi Extension of LQR
    if 1 == p % 4 or 1 == t % 4:
        return legSymb(p % t, t)
    else:
        #both p and t are === 3 mod 4
        return legSymb(p % t, t) * (-1)

#returns the numeric value of alpha mod prime given
def alphaSolver(a, prime):
    
    #extract the base and exponent
    base = a[1:a.find(')^')]
    exponent = int(a[a.find('^(') + 2:len(a) - 1])

    #from the base extract the whole number and the square root values
    whole = int(base[:base.find(' +')])
    sqaureroot = int(base[base.find('sqrt(') + 5:base.rfind(')')])

    #obtain the expression and solve using binomial expansion then simplify the result
    exp = (whole + sqrt(sqaureroot))**exponent
    eExp = expand(exp)
    sExp = simplify(eExp)

    #return only the real part mod prime
    return sExp.as_real_imag()[0] % prime

#main function of program
def main():

    parser = argparse.ArgumentParser()

    #arguments:
    #   square - the number we are trying to find the square root of
    #   prime - an odd prime
    #   verbose

    parser.add_argument('--square', action='store', help='the quadratic residue of the given prime that has a square root', required=True)
    parser.add_argument('--prime', action='store', help='an odd prime number', required=True)
    parser.add_argument('--verbose', action='store_true', help='show verbose')

    args = parser.parse_args()

    #determines of the given inputs are valid
    try:
        square = int(args.square)
        prime = int(args.prime)
    except Exception as E:
        raise E ('Only accepts Integers!')
    
    if not isPrime(prime):
        raise ValueError ('Prime number given is not a prime!')
    elif legSymb(square, prime) == -1:
        raise ValueError ('The given Square is not a quadratic residue!')

    #iterate from t = 1 and on until you reach t^2 - s is a qnr.
    t = 1
    while legSymb((t**2 - square) , prime) == 1:
        if args.verbose:
            print(f'The Legendre Symbol of ({t} / {prime}) = 1, and as such this t will be skipped')
        t += 1

    if args.verbose:
        print(f'The Legendre Symbol of ({t} / {prime}) = -1\nNow alpha will be calculated')

    alpha = f'({t} + sqrt({t**2 - square}))^({int((prime + 1) / 2)})'

    print('alpha =', alpha, '\nWhich evaluates to:', alphaSolver(alpha, prime))
    return

if __name__ == '__main__':
    try:
        main()
        
    except Exception as E:
        raise E
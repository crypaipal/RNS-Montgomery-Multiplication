import time
import sympy

def calculate_residue(number, modulus):
    return number % modulus

def rns(x, y, moduli):
    residues_x = [calculate_residue(x, modulus) for modulus in moduli]
    residues_y = [calculate_residue(y, modulus) for modulus in moduli]
    product_residues = [residues_x[i] * residues_y[i] % moduli[i] for i in range(len(moduli))]
    # for x_residue, y_residue, modulus in zip(residues_x, residues_y, moduli):
    #     product_residue = montgomeryMultiplication(x_residue, y_residue, modulus)
    #     result.append(product_residue)
    return product_residues

def extended_gcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        gcd, x, y = extended_gcd(b % a, a)
        return (gcd, y - (b // a) * x, x)


def crt(residues, moduli):
    product = 0
    prod_mod = 1

    for modulus in moduli:
        prod_mod *= modulus

    for i in range(len(residues)):
        p = prod_mod // moduli[i]
        _, inv, _ = extended_gcd(p, moduli[i])
        product += residues[i] * inv * p

    return product % prod_mod

def negative_inverse_calc(n, bit_width, base):
    def_base = 2 ** (bit_width - 1)
    neg_inverse = base - n ** (def_base - 1) % base
    return neg_inverse

def directMontgMultProduct(a, b, r, n, bit_width):
    t = a * b % r
    m = (t * negative_inverse_calc(n, bit_width, r)) % r
    s = (a * b + m * n) / r
    if s >= n:
        return s - n
    else:
        return s

def montgomeryMultiplication(x, y, moduli):
    results = []
    for modulus in moduli:
        bit_width = calc_r(modulus)
        r = 2 ** bit_width
        a = x * r % modulus
        b = y * r % modulus
        montg_product = directMontgMultProduct(a, b, r, modulus, bit_width)
        product = directMontgMultProduct(montg_product, 1, r, modulus, bit_width)
        results.append(product)
    return results

def calc_r(modul):
    for i in range(0, 32):
        r = 2 ** i
        if (r > modul) and (r / 2 <= modul):
            return i


def mod_inv(a, m):
    m0, x0, x1 = m, 0, 1
    while a > 1:
        q = a // m
        m, a = a % m, m
        x0, x1 = x1 - q * x0, x0
    return x1 + m0 if x1 < 0 else x1


def classic_MM(x, y, moduli):
    results = []
    for modulus in moduli:
        result = x * y % modulus
        results.append(result)
    #end = crt(results, moduli)
    return results

def generate_primes(num_primes):
    primes = sympy.primerange(2**15, 2**16)
    return [next(primes) for _ in range(num_primes)]

def callMenu():
    moduli = generate_primes(1000)
    print(moduli)
    while True:
        print("=====Montgomery Multiplication=====\n"
              "1 - Multiplication\n")
        try:
            choice = int(input())
        except:
            print("Input an integer number\n")
        else:
            if choice == 1:
                print("Calculator represents the x * y mod n")
                try:
                    x = int(input("Insert x: "))
                    y = int(input("Insert y: "))
                except:
                    print("x and y must be Integers!\n")
                else:
                    print('Expected (x * y mod n) FIRST result:', x * y % moduli[0])
                    start_time1 = time.perf_counter_ns()
                    result1 = rns(x, y, moduli)
                    end_time1 = time.perf_counter_ns()
                    print("The result of RNS multiplication is:", result1)
                    print("RNS Multiplication Time:", (end_time1 - start_time1), "nanoseconds\n")

                    start_time2 = time.perf_counter_ns()
                    result2 = montgomeryMultiplication(x, y, moduli)
                    end_time2 = time.perf_counter_ns()
                    print("The result of Montgomery Multiplication is:", result2)
                    print("Montgomery Multiplication Time:", (end_time2 - start_time2), "nanoseconds\n")

                    start_time3 = time.perf_counter_ns()
                    result3 = classic_MM(x, y, moduli)
                    end_time3 = time.perf_counter_ns()
                    print("The result of Classic Multiplication is:", result3)
                    print("Normal Modular Multiplication Time:", (end_time3 - start_time3), "nanoseconds\n")

                    start_time4 = time.perf_counter_ns()
                    result4 = crt(classic_MM(x, y, moduli), moduli)
                    end_time4 = time.perf_counter_ns()
                    print("The result after Chinese Remainder Theorem is: ", result4)
                    print("Chinese Remainder Theorem Time:", (end_time4 - start_time4), "nanoseconds\n")
                    holdback = input("\nPress Enter to continue:\n")
            else:
                continue

callMenu()

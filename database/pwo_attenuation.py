#!/usr/bin/env python3

import numpy
from scipy import integrate

tolerance = 1e-4


def integrand1(y, x, z, attenuation):
    return numpy.exp(-numpy.sqrt(x**2 + y**2 + z**2) * attenuation) * z / (x**2 + y**2 + z**2)**1.5


def integrand2(y, x, z, attenuation, reflectance):
    return numpy.exp(-numpy.sqrt(x**2 + y**2 + z**2) * attenuation + numpy.log(reflectance) * (x + y) / 0.0205) * z / (x**2 + y**2 + z**2)**1.5


def attenuation_fwd(depth, attenuation, reflectance):
    z = 0.18 - depth  # here depth = z_geant4 - 2.725 m

    result = 0
    n_save = 0
    for n in range(1000):
        rr = 0
        for i in range(n // 2 + 1):
            j = n - i
            x_low, x_high = 0.01025 + (i - 1) * 0.0205 if i > 0 else 0, 0.01025 + i * 0.0205
            y_low, y_high = 0.01025 + (j - 1) * 0.0205 if j > 0 else 0, 0.01025 + j * 0.0205
            r, _ = integrate.dblquad(integrand1, x_low, x_high, lambda _: y_low, lambda _: y_high, args=(z, attenuation))
            rr += r if i * 2 == n else r * 2
        rr = rr * reflectance**n
        result += rr
        if rr < result * tolerance:
            n_save = n
            break
    x_high = 0.01025 + n_save * 0.0205
    r, _ = integrate.dblquad(integrand2, 0, x_high, lambda _: x_high - _, lambda _: numpy.inf, args=(z, attenuation, reflectance))
    result += r
    r, _ = integrate.dblquad(integrand2, x_high, numpy.inf, lambda _: 0, lambda _: numpy.inf, args=(z, attenuation, reflectance))
    result += r
    return result / numpy.pi  # we only integrated in Quadrant I


def attenuation_bwd(depth, attenuation, reflectance):
    z = depth + 0.18  # here depth = z_geant4 - 2.725 m

    result = 0
    n_save = 0
    for n in range(1000):
        rr = 0
        for i in range(n // 2 + 1):
            j = n - i
            x_low, x_high = 0.01025 + (i - 1) * 0.0205 if i > 0 else 0, 0.01025 + i * 0.0205
            y_low, y_high = 0.01025 + (j - 1) * 0.0205 if j > 0 else 0, 0.01025 + j * 0.0205
            r, _ = integrate.dblquad(integrand1, x_low, x_high, lambda _: y_low, lambda _: y_high, args=(z, attenuation))
            rr += r if i * 2 == n else r * 2
        rr = rr * reflectance**(n + 1)
        result += rr
        if rr < result * tolerance:
            n_save = n
            break
    x_high = 0.01025 + n_save * 0.0205
    r, _ = integrate.dblquad(integrand2, 0, x_high, lambda _: x_high - _, lambda _: numpy.inf, args=(z, attenuation, reflectance))
    result += r * reflectance
    r, _ = integrate.dblquad(integrand2, x_high, numpy.inf, lambda _: 0, lambda _: numpy.inf, args=(z, attenuation, reflectance))
    result += r * reflectance
    return result / numpy.pi  # we only integrated in Quadrant I


if __name__ == '__main__':
    import argparse as ap

    parser = ap.ArgumentParser(prog='pwo_attenuation.py')
    parser.add_argument('inputs', nargs=2, help='penetration depth (m) and reflectance')

    args = vars(parser.parse_args())

    att_input = 1 / float(args['inputs'][0])
    ref_input = float(args['inputs'][1])

    print(att_input, ref_input)
    if ref_input > 1:
        print('bad input')
        quit()

    with open('pwo_attenuation.dat', 'w') as f:
        for depth in range(180):
            factor = attenuation_fwd(depth / 1000, att_input, ref_input) + attenuation_bwd(depth / 1000, att_input, ref_input)
            print(depth, factor)
            f.write('{}  {}\n'.format(depth, factor))

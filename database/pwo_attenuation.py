#!/usr/bin/env python3

import argparse as ap

import numpy as np
from scipy import integrate

tolerance = 1e-4


def integrand1(y, x, z, attenuation):
    r2 = x**2 + y**2 + z**2
    return np.exp(-np.sqrt(r2) * attenuation) * z / r2**1.5


def integrand2(y, x, z, attenuation, reflectance):
    r2 = x**2 + y**2 + z**2
    return (np.exp(-np.sqrt(r2) * attenuation) * reflectance**
            ((x + y) / 0.0205) * z / r2**1.5)


def attenuation_fwd(depth, attenuation, reflectance):
    z = 0.18 - depth  # here depth = z_geant4 - 2.725 m

    result = 0
    n_save = 0
    for n in range(1000):
        rr = 0
        for i in range(n // 2 + 1):
            j = n - i
            x_low = 0.01025 + (i - 1) * 0.0205 if i > 0 else 0
            x_high = 0.01025 + i * 0.0205
            y_low = 0.01025 + (j - 1) * 0.0205 if j > 0 else 0
            y_high = 0.01025 + j * 0.0205
            r, _ = integrate.dblquad(
                integrand1,
                x_low,
                x_high,
                lambda _: y_low,
                lambda _: y_high,
                args=(z, attenuation),
            )
            rr += r if i * 2 == n else r * 2
        rr = rr * reflectance**n
        result += rr
        if rr < result * tolerance:
            n_save = n
            break
    x_high = 0.01025 + n_save * 0.0205
    r, _ = integrate.dblquad(
        integrand2,
        0,
        x_high,
        lambda _: x_high - _,
        lambda _: np.inf,
        args=(z, attenuation, reflectance),
    )
    result += r
    r, _ = integrate.dblquad(
        integrand2,
        x_high,
        np.inf,
        lambda _: 0,
        lambda _: np.inf,
        args=(z, attenuation, reflectance),
    )
    result += r
    return result / np.pi  # we only integrated in Quadrant I


def attenuation_bwd(depth, attenuation, reflectance):
    z = depth + 0.18  # here depth = z_geant4 - 2.725 m

    result = 0
    n_save = 0
    for n in range(1000):
        rr = 0
        for i in range(n // 2 + 1):
            j = n - i
            x_low = 0.01025 + (i - 1) * 0.0205 if i > 0 else 0
            x_high = 0.01025 + i * 0.0205
            y_low = 0.01025 + (j - 1) * 0.0205 if j > 0 else 0
            y_high = 0.01025 + j * 0.0205
            r, _ = integrate.dblquad(
                integrand1,
                x_low,
                x_high,
                lambda _: y_low,
                lambda _: y_high,
                args=(z, attenuation),
            )
            rr += r if i * 2 == n else r * 2
        rr = rr * reflectance**(n + 1)
        result += rr
        if rr < result * tolerance:
            n_save = n
            break
    x_high = 0.01025 + n_save * 0.0205
    r, _ = integrate.dblquad(
        integrand2,
        0,
        x_high,
        lambda _: x_high - _,
        lambda _: np.inf,
        args=(z, attenuation, reflectance),
    )
    result += r * reflectance
    r, _ = integrate.dblquad(
        integrand2,
        x_high,
        np.inf,
        lambda _: 0,
        lambda _: np.inf,
        args=(z, attenuation, reflectance),
    )
    result += r * reflectance
    return result / np.pi  # we only integrated in Quadrant I


parser = ap.ArgumentParser(prog='pwo_attenuation.py')
parser.add_argument(
    'inputs', nargs=2, help='penetration depth (m) and reflectance')

args = vars(parser.parse_args())

use_ilya_method = False
if args['inputs'][0] == '-1' and args['inputs'][1] == '-1':
    use_ilya_method = True
    print('use ilya method')
else:
    att_input = 1 / float(args['inputs'][0])
    ref_input = float(args['inputs'][1])
    print(att_input, ref_input)
    if ref_input > 1 or np.isnan(att_input):
        print('bad input')
        quit()

with open('pwo_attenuation.dat', 'w') as f:
    if use_ilya_method:
        factor0 = np.exp(-1) + 1
        for dep in range(180):
            z = (dep - 90) / 180
            factor = (np.exp(-z - 0.5) + np.exp(z - 0.5)) / factor0
            print(dep, factor)
            f.write('{}  {}\n'.format(dep, factor))
    else:
        for dep in range(180):
            factor = attenuation_fwd(dep / 1000, att_input, ref_input)
            factor += attenuation_bwd(dep / 1000, att_input, ref_input)
            print(dep, factor)
            f.write('{}  {}\n'.format(dep, factor))

#!/usr/bin/env python3

import argparse as ap

import numpy as np

parser = ap.ArgumentParser(prog='hycal_module_shuffle.py')
parser.add_argument(
    'inputs', nargs=2, help='dimension uncertainties for PWO and LG')

args = vars(parser.parse_args())

err_pwo = float(args['inputs'][0]) / 1000  # um
err_lg = float(args['inputs'][1]) / 1000  # um

print(err_pwo, err_lg)

modules = np.loadtxt(
    'hycal_module.txt',
    dtype=[
        ('name', np.unicode, 10),
        ('size_x', np.double),
        ('size_y', np.double),
        ('size_z', np.double),
        ('x', np.double),
        ('y', np.double),
        ('z', np.double),
    ],
    usecols=(0, 2, 3, 4, 5, 6, 7),
)
modules['x'] = -1 * modules['x']

modules_pwo = modules[np.char.startswith(modules['name'], 'W')]
modules_lg = modules[np.char.startswith(modules['name'], 'G')]


def shuffle(module_size, n_modules, error):
    sum_size = 0
    while np.abs(sum_size - module_size * n_modules) > 1e-6:
        # notice not uniform distribution
        size = np.random.uniform(
            module_size - error,
            module_size + error,
            n_modules,
        )
        size = size * (module_size * n_modules / size.sum())
        size = size.round(4)
        sum_size = size.sum()
        if (np.any(size < module_size - error)
                or np.any(size > module_size + error)):
            sum_size = 0
    loc = np.add.accumulate(size) - size / 2
    return loc[::-1], size[::-1]


# crystal
for i in range(34):
    s = ((modules_pwo['y'] > -352.75 + (i * 20.75)) &
         (modules_pwo['y'] < -352.75 + ((i + 1) * 20.75)))
    x, size_x = shuffle(20.77, 34, err_pwo)
    if len(modules_pwo[s]) == 34:
        modules_pwo['size_x'][s] = size_x
        modules_pwo['x'][s] = x - 353.09
    else:
        modules_pwo['size_x'][s][:16] = size_x[:16]
        modules_pwo['size_x'][s][16:] = size_x[18:]
        modules_pwo['x'][s][:16] = x[:16] - 353.09
        modules_pwo['x'][s][16:] = x[18:] - 353.09

for i in range(34):
    s = ((modules_pwo['x'] > -353.09 + (i * 20.77)) &
         (modules_pwo['x'] < -353.09 + ((i + 1) * 20.77)))
    y, size_y = shuffle(20.75, 34, err_pwo)
    if len(modules_pwo[s]) == 34:
        modules_pwo['size_y'][s] = size_y
        modules_pwo['y'][s] = y - 352.75
    else:
        modules_pwo['size_y'][s][:16] = size_y[:16]
        modules_pwo['size_y'][s][16:] = size_y[18:]
        modules_pwo['y'][s][:16] = y[:16] - 352.75
        modules_pwo['y'][s][16:] = y[18:] - 352.75

# lead glass upper-left
for i in range(6):
    s = (modules_lg['x'] > -353.09) & (modules_lg['y'] > 352.75)
    s = s & ((modules_lg['y'] > 352.75 + (i * 38.15)) &
             (modules_lg['y'] < 352.75 + ((i + 1) * 38.15)))
    x, size_x = shuffle(38.15, 24, err_lg)
    modules_lg['size_x'][s] = size_x
    modules_lg['x'][s] = x - 353.09

for i in range(24):
    s = (modules_lg['x'] > -353.09) & (modules_lg['y'] > 352.75)
    s = s & ((modules_lg['x'] > -353.09 + (i * 38.15)) &
             (modules_lg['x'] < -353.09 + ((i + 1) * 38.15)))
    y, size_y = shuffle(38.15, 6, err_lg)
    modules_lg['size_y'][s] = size_y
    modules_lg['y'][s] = y - (-352.75)

# lead glass upper-right
for i in range(24):
    s = (modules_lg['x'] < -353.09) & (modules_lg['y'] > -352.75)
    s = s & ((modules_lg['y'] > -352.75 + (i * 38.15)) &
             (modules_lg['y'] < -352.75 + ((i + 1) * 38.15)))
    x, size_x = shuffle(38.15, 6, err_lg)
    modules_lg['size_x'][s] = size_x
    modules_lg['x'][s] = x - 581.99

for i in range(6):
    s = (modules_lg['x'] < -353.09) & (modules_lg['y'] > -352.75)
    s = s & ((modules_lg['x'] > -581.99 + (i * 38.15)) &
             (modules_lg['x'] < -581.99 + ((i + 1) * 38.15)))
    y, size_y = shuffle(38.15, 24, err_lg)
    modules_lg['size_y'][s] = size_y
    modules_lg['y'][s] = y - 352.75

# lead glass lower-left
for i in range(24):
    s = (modules_lg['x'] > 353.09) & (modules_lg['y'] < 352.75)
    s = s & ((modules_lg['y'] > -562.85 + (i * 38.15)) &
             (modules_lg['y'] < -562.85 + ((i + 1) * 38.15)))
    x, size_x = shuffle(38.15, 6, err_lg)
    modules_lg['size_x'][s] = size_x
    modules_lg['x'][s] = x - (-353.09)

for i in range(6):
    s = (modules_lg['x'] > 353.09) & (modules_lg['y'] < 352.75)
    s = s & ((modules_lg['x'] > 353.09 + (i * 38.15)) &
             (modules_lg['x'] < 353.09 + ((i + 1) * 38.15)))
    y, size_y = shuffle(38.15, 24, err_lg)
    modules_lg['size_y'][s] = size_y
    modules_lg['y'][s] = y - 562.85

# lead glass lower-right
for i in range(6):
    s = (modules_lg['x'] < 353.09) & (modules_lg['y'] < -352.75)
    s = s & ((modules_lg['y'] > -581.65 + (i * 38.15)) &
             (modules_lg['y'] < -581.65 + ((i + 1) * 38.15)))
    x, size_x = shuffle(38.15, 24, err_lg)
    modules_lg['size_x'][s] = size_x
    modules_lg['x'][s] = x - 562.51

for i in range(24):
    s = (modules_lg['x'] < 353.09) & (modules_lg['y'] < -352.75)
    s = s & ((modules_lg['x'] > -562.51 + (i * 38.15)) &
             (modules_lg['x'] < -562.51 + ((i + 1) * 38.15)))
    y, size_y = shuffle(38.15, 6, err_lg)
    modules_lg['size_y'][s] = size_y
    modules_lg['y'][s] = y - 581.65

modules = np.hstack((modules_lg, modules_pwo))

np.savetxt(
    'hycal_module_shuffled.dat',
    np.column_stack([
        modules['size_x'],
        modules['size_y'],
        modules['x'],
        modules['y'],
    ]),
    fmt=['%7.4f', '%7.4f', '%10.5f', '%10.5f'],
    delimiter='   ',
)

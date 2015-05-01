#! /usr/bin/env python3

import sys

def parse(filename):
    results = {}
    for line in open(filename):
        s = line.split()
        seed = int(s[0])
        np = int(s[1])
        n = int(s[2])
        score = float(s[3])
        results[seed] = {'seed': seed, 'np': np, 'n': n, 'score': score}
    return results 

def make_line(x, y):
    return '{:4d} ({:4d} {:3d}) {:>9.1f} {:>9.1f} {:>7.3f}'.format(x['seed'], x['np'], x['n'], x['score'], y['score'], x['score'] / y['score'])

np_size = sys.argv[3]
a, b = sys.argv[1] + '_' + np_size, sys.argv[2] + '_' + np_size
p, q = parse(a), parse(b)

seeds = list(set(p) & set(q))
seeds.sort(key=lambda seed: p[seed]['score'] / q[seed]['score'])

sum_ratio = 0
for seed in seeds:
    x, y = p[seed], q[seed]

    ratio = x['score'] / y['score']
    sum_ratio += ratio

    print(make_line(x, y))

total_ratio = sum_ratio / len(seeds)
print('total_ratio: {}'.format(total_ratio))

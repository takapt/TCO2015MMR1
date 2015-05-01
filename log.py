#! /usr/bin/env python3

import os
import shutil
import sys
import random
import time

HOME = os.environ['HOME']
MM = os.path.join(HOME, 'mmr1')
exe = 'a.out' # if len(sys.argv) < 2 else sys.argv[1]

output_filename = sys.argv[1]

exe_path = os.path.join(MM, exe)
copied_exe_path = os.path.join(MM, 'copied_' + str(random.randint(0, 10**5)))
shutil.copy(exe_path, copied_exe_path)

def make_command(seed):
    command = "java -jar tester.jar -exec '{}' -novis -seed {}".format(copied_exe_path, seed)
    return command

def get_score(seed):
    c = make_command(seed)

    start = time.time()
    try:
        output = os.popen(c).read()
    except:
        raise "exit"

    score = float(output.split()[-1])
    np = int(output.split()[2])
    n = int(output.split()[5])
    exe_time = time.time() - start

    return {'seed': seed, 'score': score, 'time': exe_time, 'np': np, 'n': n}

def make_line(result):
    return '{:4d} {:4d} {:3d} {:>9.1f}'.format(result['seed'], result['np'], result['n'], result['score'])

def single(seeds):
    for seed in seeds:
        result = get_score(seed)
        print('{:4d} {:3.3f} {:.3f}'.format(seed, result['score'], result['time']))
        sys.stdout.flush()

def write_file(filename, results):
    with open(filename, 'w') as f:
        for result in results:
            f.write(make_line(result) + '\n')

def multi(seeds):
    from multiprocessing import Pool
    pool = Pool()
    results = pool.map(get_score, seeds)

    small_results = []
    med_results = []
    large_results = []
    for result in results:
        if result['np'] < 100:
            small_results.append(result)
        elif result['np'] < 500:
            med_results.append(result)
        else:
            large_results.append(result)

    output_list = small_results + med_results + large_results
    for result in output_list:
        seed = result['seed']
        print(make_line(result))
        sys.stdout.flush()

    write_file(output_filename + '_s', small_results)
    write_file(output_filename + '_m', med_results)
    write_file(output_filename + '_l', large_results)

try:
#     single(range(1, 2))
    multi(range(1, 100))
finally:
    os.remove(copied_exe_path)

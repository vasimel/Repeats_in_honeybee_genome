import argparse
import os.path
from multiprocessing import Pool
from multiprocessing import cpu_count
from time import perf_counter

import dask.dataframe as dd
import pandas as pd
import tqdm

parser = argparse.ArgumentParser(description='Calculate coverage for all regions corresponding to a repeat family')
parser.add_argument('--sum_depth', default='chip_exclude/chip_exclude.txt.depth',required=True)
parser.add_argument('--trim_repeat', default='trim_honeybee.fna.out')
args = parser.parse_args()

path_to_sum_depth = args.sum_depth
path_to_trim = args.trim_repeat

trim = pd.read_csv(path_to_trim)
uniq = trim.matching_repeat.unique()


def process(nxy):
    result = sum_depth[(sum_depth[0] == nxy[0]) & (sum_depth[1].between(nxy[1], nxy[2]))][2:].drop(
        columns=[0, 1]).median()
    print(result)
    return


if os.path.exists(path_to_sum_depth + 'tmp.calc.out'):
    os.remove(path_to_sum_depth + 'tmp.calc.out')

if __name__ == '__main__':
    sum_depth = dd.read_csv(path_to_sum_depth).compute()
    start = perf_counter()
    pool = Pool(cpu_count() - 4)
    x = []
    for family in tqdm.tqdm(uniq):
        print(family)
        repeat_family_frame = trim[trim['matching_repeat'] == family]
        merged = repeat_family_frame[['query_sequence', 'pos_begin', 'pos_end']].values.tolist()
        r = pool.map(process, merged)
        x.append([family, r])
        with open(path_to_sum_depth + 'tmp.calc.out', 'a') as file:
            file.write(family + ',' + ','.join(map(str, r)) + '\n')
    pool.close()
    pool.join()
    print(perf_counter() - start)

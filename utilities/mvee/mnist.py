import pandas as pd
from sklearn.datasets import fetch_openml
from skimage.transform import resize
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import time
import argparse
import pathlib

from mvee import mvee2
from mvee import kurtosis



def read_data(col_select):
    # read in the data
    # this is somewhat slow, so we only want to do it the first time
    data_name = 'mnist_784'  # TODO
    dat = fetch_openml(data_name, as_frame=True)
    df = dat['frame']

    df_sub = df[df['class'] == df['class'].cat.categories[col_select]]
    X = df_sub.iloc[:, :-1].to_numpy()
    X = X.T
    return X


def load_cached_data(imsize, col_select):
    '''
    Read the data from a file, if available. Otherwise, read from original
    MNIST.
    '''
    fname = 'data/mnist/%d/mnist_%d.npy' % (imsize, col_select)
    try:
        X = np.load(fname)
    except Exception as e:
        print('Reading data...')
        X_big = read_data(col_select)
        print('Data read')
        X = np.zeros((imsize**2, X_big.shape[1]))
        for i in range(X.shape[1]):
            im = X_big[:, i].reshape(28, 28)
            X[:, i] = resize(im, (imsize, imsize)).flatten()
        np.save(fname, X)
    return X


def cache_data(col_select):
    fname_base = 'data/mnist/%d/mnist_%d.npy'
    X_big = read_data(col_select)
    for imsize in range(4, 20):
        fname = fname_base % (imsize, col_select)
        try:
            X = np.load(fname)
        except Exception as e:
            X = np.zeros((imsize**2, X_big.shape[1]))
            for i in range(X.shape[1]):
                im = X_big[:, i].reshape(28, 28)
                mid = 28//2
                start = mid - imsize//2
                end = start + imsize
                im = im[start:end, start:end]
                X[:, i] = im.flatten()
            np.save(fname, X)


parser = argparse.ArgumentParser(description="Run with MNIST")
parser.add_argument("method", help="{todd|newton}")
args = parser.parse_args()

if args.method not in ['todd', 'newton']:
    raise Exception('Run with argument %s. Acceptable options are '
                    '"todd" and "newton"' % args.method)

epsilon = 1e-6
all_sizes = [7]
do_todd = (args.method == 'todd')

if do_todd:
    method = 'todd'
else:
    method = 'newton'

# For writing data files
#for col in range(10):
    #cache_data(col)

kurs = np.zeros((10, len(all_sizes)))
times = np.zeros((10, len(all_sizes)))
iters = np.zeros((10, len(all_sizes)))
cores = np.zeros((10, len(all_sizes)))
for digit in range(10):
    print('digit: %d' % digit)
    for imsize in all_sizes:
        X = load_cached_data(imsize, digit)
        print('m: ', X.shape[1])
        kurs[digit, all_sizes.index(imsize)] = kurtosis(X, do_log=False)
        try:
            t1 = time.time()
            ret = mvee2(X, initialize='qr', epsilon=epsilon, method=method,
                        verbose=True, max_iter=100000, full_output=True,
                        upproject=True, track_count=True)
            t2 = time.time()
            times[digit, all_sizes.index(imsize)] = t2 - t1
            iters[digit, all_sizes.index(imsize)] = ret['iter_count']
            cores[digit, all_sizes.index(imsize)] = (ret['u'] > 1e-12).sum()
            print('total time: ', t2 - t1)
            print('iterations: ', ret['iter_count'])
            print('')
        except:
            times[digit, all_sizes.index(imsize):] = np.inf
            iters[digit, all_sizes.index(imsize):] = np.inf
            cores[digit, all_sizes.index(imsize):] = np.inf
            break

#firstinf = np.where(times == np.inf)[1].min()
firstinf = kurs.shape[1]
skurs = kurs[:, :firstinf]
stimes = times[:, :firstinf]
siters = iters[:, :firstinf]

base_dir = 'outputs/mnist/'
if do_todd:
    full_dir = base_dir + 'coord-ascent/'
else:
    full_dir = base_dir + 'newton/'
pathlib.Path(full_dir).mkdir(parents=True, exist_ok=True)
np.save(full_dir + 'kurtosis', kurs)
np.save(full_dir + 'times', times)
np.save(full_dir + 'iters', iters)
np.save(full_dir + 'cores', cores)

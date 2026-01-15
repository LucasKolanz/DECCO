'''
Run MVEE on the document dump data set. The document word counts are stored in
a format similar to COO. Choose n words as our number of dimensions.
The number of documents is m. Each column in a single matrix is the word
counts for a single document for some size-n subset of the
words.

generate_matrix and read_data can be used for smaller problems where the entire
data set fits into memory, but they may take up a large amount of disk space
when writing the matrix.

For moderately large problems, generate_submatrix should be used because it
avoids turning the entire COO matrix into one large dense matrix and instead
creates only as many rows as needed.

If matrices were generated with generate_submatrix, they can be read all at
once with read_matrices (convenient for small problems), or one at a time
with read_matrix (useful when the problem is large).
'''
import pandas as pd
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mvee import mvee2
from mvee import kurtosis
import time
import random
import scipy.linalg as sla
import argparse
import pathlib



def generate_matrix(fname):
    '''
    Read the COO-like file into a matrix and save that matrix to disk to
    avoid needing to do this every time.

    Note: This does not work for larger data sets.
    '''
    with open('data/doc-dumps/docword.%s.txt' % fname, 'r') as f:
        m = int(f.readline())
        n = int(f.readline())
        mat = np.zeros((n, m))
        f.readline()
        for line in f:
            j, i, d = [int(l) for l in line.split()]
            i -= 1
            j -= 1
            mat[i, j] = d
    np.save('data/doc-dumps/matrices/%s' % fname, mat)
    return mat


def generate_submatrix(fname, n, n_runs):
    '''
    Generate a matrix using only the specified rows.
    '''
    print('Generating matrices...')
    mult_factor = 4  # needs to be larger to prevent rank deficiency
    with open('data/doc-dumps/docword.%s.txt' % fname, 'r') as f:
        m = int(f.readline())

        rows = []
        mats = []
        for i in range(n_runs):
            # TODO grab more rows than needed
            rows.append(random.sample(range(m), mult_factor*n))
            mats.append(np.zeros((mult_factor*n, m)))

        f.readline()  # ignore size of dictionary
        total_rows = int(f.readline())
        count = 0
        for line in f:
            count += 1
            if count % 100000 == 0:
                print('Line %d of %d' % (count, total_rows))
            j, i, d = [int(l) for l in line.split()]
            i -= 1
            j -= 1
            for ind in range(n_runs):
                if i in rows[ind]:
                    mats[ind][rows[ind].index(i), j] = d
    
    # Remove extra rows
    for i in range(len(mats)):
        mat = mats[i]
        _, _, P = sla.qr(mat.T, mode='economic', pivoting=True, check_finite=False)
        mats[i] = mat[P[:(n + 1)]]  # (n + 1) because of upproject
        

    base_mat_dir = 'data/doc-dumps/matrices/%s/' % fname
    for i in range(n_runs):
        np.save(base_mat_dir + 'mat_%d' % i, mats[i])

    if la.matrix_rank(mats[0]) < n:
        raise Exception
    print('Finished generating matrices.')
    return m, mats


def read_data(fname, n):
    '''
    Note: this does not work for the larger data sets because they cannot
    easily be stored as a dense matrix.
    '''
    print('Reading data...')
    try:
        mats = []
        for val in ['min', 'med', 'max']:
            mats.append(np.load('data/doc-dumps/matrices/ordered/%s/%s.npy' %
                        (fname, val)))
        m = mats[0].shape[1]
    except:
        print('Generating matrix...')
        m, mats = generate_submatrix(fname, n)
        print('Finished generating matrix.')
    print('Data read.')
    return m, mats


def read_matrix(fname, n, i):
    print('Reading matrix...')

    base_mat_dir = 'data/doc-dumps/matrices/ordered/%s/' % fname
    fname = base_mat_dir + 'mat_%d.npy' % i
    mat = np.load(fname)

    print('Finished reading matrix.')
    return mat.shape[1], mat



parser = argparse.ArgumentParser(description="Generate doc-dump results")
parser.add_argument("method", help="{todd|newton}")
parser.add_argument("fname", help="{enron|kos|nips|nytimes}")
args = parser.parse_args()

if args.method not in ['todd', 'newton']:
    raise Exception('Run with argument %s. Acceptable options are '
                    '"todd" and "newton"' % args.method)

if args.fname not in ['enron', 'kos', 'nips', 'nytimes']:
    raise Exception('%s is not a valid data set')

do_todd = (args.method == 'todd')
fname = args.fname
epsilon = 1e-6
n = 50
n_runs = 3

if do_todd:
    added_dir = 'coord-ascent/'
    method = 'todd'
    max_iter = 10000
else:
    added_dir = ''
    method = 'newton'
    max_iter = 200

base_dir = './outputs/doc-dumps/%s%s/' % (added_dir, fname)
pathlib.Path(base_dir).mkdir(parents=True, exist_ok=True)

m, mats = read_data(fname, n)

try:
    times = np.load(base_dir + 'times.npy')
    cores = np.load(base_dir + 'cores.npy')
    iters = np.load(base_dir + 'iters.npy')
    kurs = np.load(base_dir + 'kurtosis.npy')
    start = np.argmin(times > 0)
except:
    times = np.zeros(n_runs)
    kurs = np.zeros(n_runs)
    cores = np.zeros(n_runs)
    iters = np.zeros(n_runs)
    start = 0
for i in range(start, n_runs):
    print('Starting iteration %d of %d' % (i + 1, n_runs))
    X = mats[i]
    kurs[i] = kurtosis(X, do_log=False)

    try:
        t1 = time.time()
        ret = mvee2(X, initialize='qr', epsilon=epsilon, method=method,
                    verbose=True, max_iter=max_iter, full_output=True,
                    upproject=True, track_count=True)
        t2 = time.time()
        times[i] = t2 - t1
        iters[i] = ret['iter_count']
        cores[i] = (ret['u'] > 1e-12).sum()
        print('total time: ', t2 - t1)
        print('iterations: ', ret['iter_count'])
        print('')
    except Exception as e:
        times[i] = np.inf
        iters[i] = np.inf
        cores[i] = np.inf

    # Data must be saved every iteration because the program tends to get
    # killed when running on larger data
    np.save(base_dir + 'times', times)
    np.save(base_dir + 'cores', cores)
    np.save(base_dir + 'iters', iters)
    np.save(base_dir + 'kurtosis', kurs)

if 0:
    plt.figure()
    plt.semilogx(kurs, cores, 'o')
    #plt.semilogx(kurs, iters, 'o')
    plt.title('Document data (%s)\nm = %d; n = %d' % (fname, m, n))
    plt.xlabel('kurtosis')
    #plt.ylabel('iterations')
    plt.ylabel('core set size')
    #plt.savefig('/home/nate/show-mike-images/doc-dumps/enron')
    #plt.close('all')
    plt.show()

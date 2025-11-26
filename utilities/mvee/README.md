This software computes the minimum-volume enclosing ellipsoid (MVEE) of a
data set.
The code for MVEE computation is in the file `mvee.py`.
We suggest using the function `mvee2` with the `hybrid` method for good
performance across a broad array of data sets.
Details about the code and its use are provided in the paper referenced below.
For a brief overview of the code, see `CODE-STRUCTURE.txt`.

Other than `mvee.py`,
the remainder of the code exists to facilitate reproduction of the results from
the document cited below.
For advice on running the experiments in that document, see `RUNNING-CODE.txt`.

# License

The software is released open-source under the MIT license.

The data sets are taken from the UC Irvine Machine Learning Repository and
were made available under the
Creative Commons Attribution 4.0 International license (CC BY 4.0).

# Citing

If you make use of this software in a publication,
please cite:

```
Bowman, N.; Heath, M. T., "Computing Minimum-Volume Enclosing Ellipsoids,"
Mathematical Programming Computation, to appear
```

If you use the data sets (i.e., files under `data/`) in a publication,
please cite

```
Dua, D. and Graff, C. (2019). UCI Machine Learning Repository [http://archive.ics.uci.edu/ml]. Irvine, CA: University of California, School of Information and Computer Science.
```

import h5py
f = h5py.File('example.h5','r')
for item in f.keys():
    for it in f[item]:
        print(it)
        print (item + ":", f['/'+item][it][:])
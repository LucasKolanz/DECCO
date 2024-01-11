import h5py

file = '/home/lucas/Desktop/SpaceLab_data/jobs/error2Test2/N_10/T_10/7_data.h5'

# f = h5py.File('/home/lucas/Desktop/SpaceLab_data/restartTest1/N_5/T_3/2_data.h5','r')
f = h5py.File(file,'r')

data = f['/simData'][-2000:-1000]
print(data)
print(len(data))


# for item in f.keys():
#     print(item)
#     for it in f[item]:
#         print(it)
#         print (item + ":", f['/'+item])
        # print (item + ":", f['/'+item][it][:])





# Replace 'file_path' with your HDF5 file path and 'dataset_name' with your dataset name
# file_path = '/home/lucas/Desktop/SpaceLab_data/jobs/restartTest1/N_5/T_3/3_data.h5'
# # file_path = '/home/lucas/Desktop/SpaceLab_data/jobs/restartTest1/N_5/T_3/2_data.h5'
# dataset_name = 'writes'

# with h5py.File(file_path, 'r') as file:
#     dataset = file[dataset_name]
#     # Assuming the value you want to read is at a specific index, for example [0,0] for a 2D dataset
#     scalar_value = dataset[()]

# print(scalar_value)
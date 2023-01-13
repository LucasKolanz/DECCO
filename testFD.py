import utils as u
# from scipy.spatial.transform import Rotation as R
import numpy as np

# data_folder = '/mnt/be2a0173-321f-4b9d-b05a-addba547276f/kolanzl/SpaceLab/jobs/comparison_test1/N_100/T_1000/'
# o3dv = u.o3doctree(data_folder,overwrite_data=True)
# o3dv.make_tree()


def rotation_matrix(v1, v2):
    """
    Returns the rotation matrix between two vectors v1 and v2.
    Both v1 and v2 must be numpy arrays with the same shape.

    :param v1: First vector
    :param v2: Second vector
    :return: Rotation matrix
    """
    v1 = np.array(v1)
    v2 = np.array(v2)
    if v1.shape != v2.shape:
        raise ValueError("Both vectors must have the same shape.")
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    v = np.cross(v1, v2)
    s = np.linalg.norm(v)
    c = np.dot(v1, v2)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(v1.shape[0]) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
    return rotation_matrix



v1 = [1, 0, 0]
v2 = [0, 1, 0]
expected = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
print(expected)
print(rotation_matrix(v1, v2))
# np.testing.assert_almost_equal(rotation_matrix(v1, v2), expected)

# v1 = [1, 2, 3]
# v2 = [2, 3, 4]
# expected = [[-0.04641016, -0.89442719,  0.4472136 ],\
#  [ 0.9916198 , -0.09053575, -0.09053575],\
#  [ 0.11952136,  0.4472136 ,  0.88191711]]
# np.testing.assert_almost_equal(rotation_matrix(v1, v2), expected, decimal=6)
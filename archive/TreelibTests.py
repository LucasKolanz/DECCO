# from treelib import Node, Tree
# import uuid
import numpy as np


# tree = Tree()

stri = '001'
strj = '101'

# tree.create_node(uuid.uuid1(),data=stri)


class node(object):
	"""docstring for node"""
	def __init__(self,parent,id,data):
		super(node, self).__init__()
		self.parent = parent
		self.data = data
		self.id = uuid.uuid1()
		self.kids = []


class tree(object):
	"""docstring for tree"""
	def __init__(self, root):
		super(tree, self).__init__()
		self.root = root
		self.nnodes = 0 #root is node 0
	
	# def add_node(parent,id,data):



x, y, z = np.indices((8, 8, 8))

# draw cuboids in the top left and bottom right corners, and a link between
# them
cube1 = (x < 3) & (y < 3) & (z < 3)
cube2 = (x >= 5) & (y >= 5) & (z >= 5)
link = abs(x - y) + abs(y - z) + abs(z - x) <= 2

# combine the objects into a single boolean array
voxelarray = cube1 | cube2 | link

print(x)

exit(0)

result = 0

if result == 1:
	print('bad')
elif result == 0:
	print("bad1")
else:
	print('good')
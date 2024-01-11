import ball_group as bg
import matplotlib.pyplot as plt

def sim_looper(O):
	for Step in range(O.steps):
		writeStep = False
		if Step % O.skip == 0:
			print("step: {}".format(Step))
			writeStep = True
		O.sim_one_step(writeStep)
		# plt.show()
		# exit(0)

def main():
	O = bg.Ball_group(4)
	# O.sim_one_step(init=True)
	sim_looper(O)

	# for i in range(O.num_balls):
	# 	O.add_projectile()
	# 	sim_looper()

if __name__ == '__main__':
	main()
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Ball_group:
	"""docstring for Ball_group"""
	def __init__(self, num_particles=4):
		super(Ball_group, self).__init__()

		#plot stuff
		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(111, projection='3d')
		self.ax.set_box_aspect([1,1,1])
		
		#sim stuff
		self.num_particles = num_particles

		self.mom = [0.0,0.0,0.0] 
		self.ang_mom = [0.0,0.0,0.0]
		self.PE = 0.0
		self.KE = 0.0

		# NEED TO SET
		self.dt = 6.45497e-11
		self.kout = 11.3097
		self.kin = 28.2743
		self.h_min = 1e-6
		self.Ha = 4.7e-12
		self.u_s = 0.1
		self.u_r = 1e-5
		self.steps = 7745966
		self.skip = 154919

		self.distances = np.zeros((int((num_particles * num_particles / 2) - (num_particles / 2))))
		self.pos = np.zeros((num_particles,3))
		self.vel = np.zeros((num_particles,3))
		self.velh = np.zeros((num_particles,3))
		self.acc = np.zeros((num_particles,3))
		self.w = np.zeros((num_particles,3))
		self.wh = np.zeros((num_particles,3))
		self.aacc = np.zeros((num_particles,3))

		# self.R = np.zeros((num_particles))
		# self.m = np.zeros((num_particles))
		# self.moi = np.zeros((num_particles))

		self.R = np.array([1e-5,1e-5,5e-6,5e-6])
		self.m = np.array([9.42478e-15,9.42478e-15,1.1781e-15,1.1781e-15])
		self.moi = np.array([3.76991e-25,3.76991e-25,1.1781e-26,1.1781e-26])

		self.pos[0] = np.array([-7.05566e-06,9.81979e-06,1.17522e-06])
		self.pos[1] = np.array([1.36658e-06,-7.71631e-06,-1.76359e-06])
		self.pos[2] = np.array([1.24614e-06,1.29462e-06,9.97946e-06])
		self.pos[3] = np.array([4.42665e-05,-1.81225e-05,-5.27251e-06])

		self.vel[0] = np.array([0.012564,-0.0185357,0.00228819])
		self.vel[1] = np.array([0.0204142,0.00294446,-0.00275746])
		self.vel[2] = np.array([-0.00130849,-0.030628,-0.0244852])
		self.vel[3] = np.array([-0.262517,0.155358,0.0282393])

		self.plot_spheres()



	def sim_one_step(self,write_step=False,init=False):
		if not init:
			for Ball in range(self.num_particles):
				self.velh[Ball] = self.vel[Ball] + 0.5*self.acc[Ball]*self.dt
				self.wh[Ball] = self.w[Ball] + 0.5*self.aacc[Ball]*self.dt
				self.pos[Ball] += self.velh[Ball]*self.dt
				self.acc[Ball] = np.array([0,0,0])
				self.aacc[Ball] = np.array([0,0,0])

		for A in range(1,self.num_particles):
			for B in range(0,A):
				sumRaRb = self.R[A] + self.R[B]
				rVecab = self.pos[B] - self.pos[A]
				rVecba = -rVecab
				dist = la.norm(rVecab)

				overlap = sumRaRb-dist

				totalForceOnA = np.array([0.0,0.0,0.0])

				e = int((A * (A-1) * 0.5) + B)
				oldDist = self.distances[e]

				if (overlap > 0):
					if dist >= oldDist:
						k = self.kout
					else:
						k = self.kin

					h = self.h_min
					Ra = self.R[A]
					Rb = self.R[B]
					h2 = h*h
					twoRah = 2*Ra*h
					twoRbh = 2*Rb*h
					vdwForceOnA = self.Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb * \
							 ((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) * \
											   (h2 + twoRah + twoRbh + 4 * Ra * Rb) * \
											   (h2 + twoRah + twoRbh + 4 * Ra * Rb))) * \
							 la.norm(rVecab);

					elasticforceOnA = -k * overlap * 0.5 * (rVecab/dist)

					slideForceOnA = np.zeros((3))
					rollForceOnA = np.zeros((3))
					# torqueA = np.zeros((3))
					# torqueB = np.zeros((3))

					# Shared terms;
					elastic_force_A_mag = la.norm(elasticforceOnA)
					r_a = rVecab * Ra/sumRaRb
					r_b = rVecba * Rb/sumRaRb
					w_diff = self.w[A] - self.w[B]

					# Sliding fric terms
					d_vel = self.vel[B] - self.vel[A]
					frame_A_vel_B = d_vel - np.dot(d_vel,rVecab) * (rVecab/(dist*dist)) - \
									np.cross(self.w[A],r_a) - np.cross(self.w[B],r_a)

					# Compute sliding fric force
					rel_vel_mag = la.norm(frame_A_vel_B)
					if rel_vel_mag > 1e-13:
						slideForceOnA = self.u_s * elastic_force_A_mag * (frame_A_vel_B/rel_vel_mag)

					#Compute rolling fric force
					w_diff_mag = la.norm(w_diff)
					if w_diff_mag > 1e-13:
						rollForceOnA = -self.u_r * elastic_force_A_mag * \
									(np.cross(w_diff,r_a)/la.norm(np.cross(w_diff,r_a)))

					#Total forces on a:
					# ?????????????????????? SHould this add on rolling force ????????????????
					totalForceOnA = elasticforceOnA + slideForceOnA + vdwForceOnA

					#Total torque on a and b:
					torqueA = np.cross(r_a, slideForceOnA + rollForceOnA)
					torqueB = np.cross(r_b, -slideForceOnA + rollForceOnA)

					self.aacc[A] += torqueA / self.moi[A]
					self.aacc[B] += torqueB / self.moi[B]


					if write_step:
						self.PE += 0.5 * k * overlap * overlap

						diffRaRb = self.R[A] - self.R[B] 
						z = sumRaRb + h
						two_RaRb = 2 * self.R[A] * self.R[B]
						denom_sum = z * z - (sumRaRb*sumRaRb)
						denom_diff = z * z - (diffRaRb*diffRaRb)
						U_vdw = -(self.Ha/6) * (two_RaRb / denom_sum + two_RaRb / denom_diff + np.log(denom_sum / denom_diff))

						self.PE += U_vdw

				else:
					h = np.absolute(overlap)
					if h < self.h_min:
						h = self.h_min

					Ra = self.R[A]
					Rb = self.R[B]
					h2 = h*h
					twoRah = 2 * Ra * h
					twoRbh = 2 * Rb * h
					vdwForceOnA = self.Ha / 6 * 64 * Ra * Ra * Ra * Rb * Rb * Rb * \
										((h + Ra + Rb) / ((h2 + twoRah + twoRbh) * (h2 + twoRah + twoRbh) * \
										(h2 + twoRah + twoRbh + 4 * Ra * Rb) * \
										(h2 + twoRah + twoRbh + 4 * Ra * Rb))) * \
										la.norm(rVecab);

					totalForceOnA = vdwForceOnA

					if write_step:
						diffRaRb = self.R[A] - self.R[B] 
						z = sumRaRb + h
						two_RaRb = 2 * self.R[A] * self.R[B]
						denom_sum = z * z - (sumRaRb*sumRaRb)
						denom_diff = z * z - (diffRaRb*diffRaRb)
						U_vdw = -(self.Ha/6) * (two_RaRb / denom_sum + two_RaRb / denom_diff + np.log(denom_sum / denom_diff))

						self.PE += U_vdw

				self.acc[A] += totalForceOnA / self.m[A] 
				self.acc[B] -= totalForceOnA / self.m[B] 

				self.distances[e] = dist


		#Third pass - calculate next step velocity
		for Ball in range(self.num_particles):
			self.vel[Ball] = self.velh[Ball] + 0.5*self.acc[Ball]*self.dt
			self.w[Ball] = self.wh[Ball] + 0.5*self.aacc[Ball]*self.dt

			if write_step:

				self.KE += 0.5 * self.m[Ball] * np.inner(self.vel[Ball],self.vel[Ball]) + \
								0.5 * self.moi[Ball]*np.inner(self.w[Ball],self.w[Ball])
				self.mom += self.m[Ball]*self.vel[Ball]
				self.ang_mom += self.m[Ball] * np.cross(self.pos[Ball],self.vel[Ball]) + self.moi[Ball] * self.w[Ball]

		if write_step:
			self.plot_spheres()


	def plot_spheres(self):
		
		self.ax.cla()
		# Plot the spheres
		max_axis = 0.0
		min_axis = 0.0
		for i in range(len(self.pos)):
			u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
			x = self.pos[i][0] + self.R[i] * np.cos(u) * np.sin(v)
			y = self.pos[i][1] + self.R[i] * np.sin(u) * np.sin(v)
			z = self.pos[i][2] + self.R[i] * np.cos(v)
			if np.max(x) > max_axis:
				max_axis = np.max(x)
			elif np.max(y) > max_axis:
				max_axis = np.max(y)
			elif np.max(z) > max_axis:
				max_axis = np.max(z)

			if np.min(x) < min_axis:
				min_axis = np.min(x)
			elif np.min(y) < min_axis:
				min_axis = np.min(y)
			elif np.min(z) < min_axis:
				min_axis = np.min(z)
			self.ax.plot_surface(x, y, z, color='b', alpha=0.6)

		# Set the axis labels and limits
		self.ax.set_xlabel('X')
		self.ax.set_ylabel('Y')
		self.ax.set_zlabel('Z')
		self.ax.set_xlim([min_axis, max_axis])
		self.ax.set_ylim([min_axis, max_axis])
		self.ax.set_zlim([min_axis, max_axis])

		# Show the plot
		plt.pause(0.0001)
		# plt.show()





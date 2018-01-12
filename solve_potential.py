import numpy as np 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt 
from constant import *
import random, time, sys
from skimage import io

def sign(x):	return 1 if x>=0 else -1

class Particle(object):
	def __init__(self, p, v, m):
		self.p = p
		self.v = v
		self.a = np.zeros(2)
		self.m = m
	'''
	def update(self):
		t = 1
		self.v += self.a * t
		self.p += self.v * t
		self.track.append(copy.copy(self.p))
	
	def draw(self, surface):
		pygame.draw.circle(surface, self.c, [int(self.p[0]), int(self.p[1])], 10)
		for each in self.track:
			surface.fill(self.c, (each,(2,2)))
	'''


class Grid(object):
	"""docstring for Grid"""

	def __init__(self, n):
		self.n = n
		self.particles = [Particle(np.array([random.random()*SCREENSIZE[0],random.random()*SCREENSIZE[1]]),
							np.array([0.0, 0.0]), 10**8) for i in range(n)]
		# density on the grid
		self.grid = np.zeros(SCREENSIZE)
		# acceleration on the grid
		self.ax = np.zeros(SCREENSIZE)
		self.ay = np.zeros(SCREENSIZE)
		# Green function and its fourier transform
		Green = np.array([np.array([ 1/np.sqrt((i - SCREENSIZE[0]/2)**2 + (j - SCREENSIZE[1]/2)**2) 
										for j in range(SCREENSIZE[1])])
										 for i in range(SCREENSIZE[0])])
		Green[SCREENSIZE[0]/2,SCREENSIZE[1]/2] = 1
		self.G_fft = np.fft.fft2(Green)
		


		

	def particle_to_density(self):
		"""density of grid is interpolated using cloud-in-cell method"""
		self.grid = np.zeros(SCREENSIZE)
		for each in self.particles:
			x = each.p
			m = each.m
			i, j = int(x[0]), int(x[1])
			xi, xj = x[0] - i, x[1] - j
			assert xi > 0, "xi < 0"
			assert xj > 0, "xj < 0, {0}".format(x[1])

			xi-=0.5; xj-=0.5
			self.grid[i,j] 						+= m * (1 - abs(xi)) * (1 - abs(xj))
			self.grid[(i+sign(xi))%SCREENSIZE[0],j] 	+= m * abs(xi) * (1 - abs(xj))
			self.grid[i,(j+sign(xj))%SCREENSIZE[1]]	+= m * (1 - abs(xi)) * abs(xj)
			self.grid[(i+sign(xi))%SCREENSIZE[0],(j+sign(xj))%SCREENSIZE[1]] 	+= m * abs(xi) * abs(xj)

			'''
			if xi > 0.5:
				if xj > 0.5:
					self.grid[i,j] 						+= m * (1.5 - xi) * (1.5 - xj)
					self.grid[(i+1)%SCREENSIZE[0],j] 	+= m * (xi - 0.5) * (1.5 - xj)
					self.grid[i,(j+1)%SCREENSIZE[1]]	+= m * (1.5 - xi) * (xj - 0.5)
					self.grid[(i+1)%SCREENSIZE[0],(j+1)%SCREENSIZE[1]] 	+= m * (xi - 0.5) * (xj - 0.5)

				if xj <= 0.5:
					self.grid[i,j] 		+= m * (1.5 - xi) * (0.5 + xj)
					self.grid[(i+1)%SCREENSIZE[0],j] 	+= m * (xi - 0.5) * (0.5 + xj)
					self.grid[i,j-1]	+= m * (1.5 - xi) * (0.5 - xj)
					self.grid[(i+1)%SCREENSIZE[0],j-1] 	+= m * (xi - 0.5) * (0.5 - xj)
			if xi <= 0.5:
				if xj > 0.5:
					self.grid[i,j] 		+= m * (xi + 0.5) * (1.5 - xj)
					self.grid[i-1,j] 	+= m * (0.5 - xi) * (1.5 - xj)
					self.grid[i,(j+1)%SCREENSIZE[1]]	+= m * (xi + 0.5) * (xj - 0.5)
					self.grid[i-1,(j+1)%SCREENSIZE[1]] 	+= m * (0.5 - xi) * (xj - 0.5)

				if xj <= 0.5:
					self.grid[i,j] 		+= m * (xi + 0.5) * (0.5 + xj)
					self.grid[i-1,j] 	+= m * (0.5 - xi) * (0.5 + xj)
					self.grid[i,j-1]	+= m * (xi + 0.5) * (0.5 - xj)
					self.grid[i-1,j-1] 	+= m * (0.5 - xi) * (0.5 - xj)
			'''

	def calculate_potential(self):
		""" 
			poisson equation:
			d^2/dx^2 \phi(x) = 4 \pi G \rho(x)

			Green function:
			d^2/dx^2 Green(x) = \delta(x)

			solution of poisson equation is convolution of Green function and density function:
			\phi(x) = 4 \pi G (Green * \rho)(x)

			fourier transform:
			\phi(k) = 4 \pi G Green(k) * \rho(k)
			
		"""
		density_fft = 4 * np.pi * G * np.fft.fft2(self.grid)
		self.potential = np.fft.fftshift(np.fft.ifft2(density_fft * self.G_fft)).real


	def calculate_force(self):
		"""
			4th order differentiation of potential
			f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
			F = -\Delta \phi * m

			interpolating using the same cloud-in-cell method as assigning particles to grid
		"""
		factor = [-1.0/12, 2.0/3, 0, -2.0/3, 1.0/12]
		for j in range(SCREENSIZE[1]):
			self.ax[:,j] = np.convolve(self.potential[:,j], factor, "same")
		for i in range(SCREENSIZE[0]):
			self.ay[i,:] = np.convolve(self.potential[i,:], factor, "same")
		'''
		for i in range(SCREENSIZE[0]):
			for j in range(SCREENSIZE[1]):
				self.ax[i,j] = (-self.potential[(i+2) % SCREENSIZE[0], j] + 8 * self.potential[(i+1) % SCREENSIZE[0], j]
								-8 * self.potential[i-1, j] + self.potential[i-2, j]) / 12
				self.ay[i,j] = (-self.potential[i, (j+2) % SCREENSIZE[1]] + 8 * self.potential[i, (j+1) % SCREENSIZE[1]]
								-8 * self.potential[i, j-1] + self.potential[i, j-2]) / 12
		'''
		for each in self.particles:
			x = each.p
			i, j = int(x[0]), int(x[1])
			xi, xj = x[0] - i, x[1] - j
			xi-=0.5; xj-=0.5
			inext, jnext = (i+sign(xi))%SCREENSIZE[0], (j+sign(xj))%SCREENSIZE[1]
			each.a[0] = self.ax[i,j] * (1 - abs(xi)) * (1 - abs(xj)) + \
						self.ax[inext,j] * abs(xi) * (1 - abs(xj)) + \
						self.ax[i,jnext] * (1 - abs(xi)) * abs(xj) + \
						self.ax[inext,jnext] * abs(xi) * abs(xj)
			each.a[1] = self.ay[i,j] * (1 - abs(xi)) * (1 - abs(xj)) + \
						self.ay[inext,j] * abs(xi) * (1 - abs(xj)) + \
						self.ay[i,jnext] * (1 - abs(xi)) * abs(xj) + \
						self.ay[inext,jnext] * abs(xi) * abs(xj)


			'''
			if xi > 0.5:
				if xj > 0.5:
					each.a[0] = self.ax[i,j] * (1.5 - xi) * (1.5 - xj) + \
								self.ax[(i+1)%SCREENSIZE[0],j] * (xi - 0.5) * (1.5 - xj) + \
								self.ax[i,(j+1)%SCREENSIZE[1]] * (1.5 - xi) * (xj - 0.5) + \
								self.ax[(i+1)%SCREENSIZE[0],(j+1)%SCREENSIZE[1]] * (xi - 0.5) * (xj - 0.5)
					each.a[1] = self.ay[i,j] * (1.5 - xi) * (1.5 - xj) + \
								self.ay[(i+1)%SCREENSIZE[0],j] * (xi - 0.5) * (1.5 - xj) + \
								self.ay[i,(j+1)%SCREENSIZE[1]] * (1.5 - xi) * (xj - 0.5) + \
								self.ay[(i+1)%SCREENSIZE[0],(j+1)%SCREENSIZE[1]] * (xi - 0.5) * (xj - 0.5)
				if xj <= 0.5:
					each.a[0] = self.ax[i,j] * (1.5 - xi) * (0.5 + xj) +\
								self.ax[(i+1)%SCREENSIZE[0],j]  * (xi - 0.5) * (0.5 + xj) +\
								self.ax[i,j-1] * (1.5 - xi) * (0.5 - xj) +\
								self.ax[(i+1)%SCREENSIZE[0],j-1]  * (xi - 0.5) * (0.5 - xj)
					each.a[1] = self.ay[i,j] * (1.5 - xi) * (0.5 + xj) +\
								self.ay[(i+1)%SCREENSIZE[0],j]  * (xi - 0.5) * (0.5 + xj) +\
								self.ay[i,j-1] * (1.5 - xi) * (0.5 - xj) +\
								self.ay[(i+1)%SCREENSIZE[0],j-1]  * (xi - 0.5) * (0.5 - xj)

			if xi <= 0.5:
				if xj > 0.5:
					each.a[0] = self.ax[i,j] * (xi + 0.5) * (1.5 - xj) +\
								self.ax[i-1,j] * (0.5 - xi) * (1.5 - xj) +\
								self.ax[i,(j+1)%SCREENSIZE[1]] * (xi + 0.5) * (xj - 0.5) +\
								self.ax[i-1,(j+1)%SCREENSIZE[1]] * (0.5 - xi) * (xj - 0.5)
					each.a[1] = self.ay[i,j] * (xi + 0.5) * (1.5 - xj) +\
								self.ay[i-1,j] * (0.5 - xi) * (1.5 - xj) +\
								self.ay[i,(j+1)%SCREENSIZE[1]] * (xi + 0.5) * (xj - 0.5) +\
								self.ay[i-1,(j+1)%SCREENSIZE[1]] * (0.5 - xi) * (xj - 0.5)

				if xj <= 0.5:
					each.a[0] = self.ax[i,j] * (xi + 0.5) * (0.5 + xj) +\
								self.ax[i-1,j] * (0.5 - xi) * (0.5 + xj) +\
								self.ax[i,j-1] * (xi + 0.5) * (0.5 - xj) +\
								self.ax[i-1,j-1] * (0.5 - xi) * (0.5 - xj)
					each.a[1] = self.ay[i,j] * (xi + 0.5) * (0.5 + xj) +\
								self.ay[i-1,j] * (0.5 - xi) * (0.5 + xj) +\
								self.ay[i,j-1] * (xi + 0.5) * (0.5 - xj) +\
								self.ay[i-1,j-1] * (0.5 - xi) * (0.5 - xj)
				'''

	def update(self):
		t = 1
		for each in self.particles:
			each.p += each.v * t
			each.p[0] %= SCREENSIZE[0]
			each.p[1] %= SCREENSIZE[1]
			each.v += each.a * t

	def draw(self):
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		X = np.arange(0, SCREENSIZE[1])
		Y = np.arange(0, SCREENSIZE[0])
		X, Y = np.meshgrid(X, Y)
		surf = ax.plot_wireframe(X, Y, self.potential, rstride = 20, cstride = 20, cmap = cm.coolwarm, linewidth=1)
		#ax.set_zlim(320, 350)
		#fig.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()

if __name__ == "__main__":
	start_time = time.time()
	a = Grid(10**5)
	#print sys.getsizeof(a.particles)/1024/1024
	print "init time: {0}".format(time.time() - start_time)
	for i in xrange(100):
		a.particle_to_density()
		print "particle_to_density time: {0}".format(time.time() - start_time)
		plt.figure()
		plt.imshow(a.grid.T/10**8, cmap = 'magma', vmin = 0, vmax = 1)
		plt.axis('off')
		plt.savefig('{0}.png'.format(i))
		plt.close()
		print "plot time: {0}".format(time.time() - start_time)
		a.calculate_potential()
		print "calculate_potential time: {0}".format(time.time() - start_time)
		a.calculate_force()
		print "calculate_force time: {0}".format(time.time() - start_time)
		a.update()
	#a.draw()


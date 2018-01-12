import pygame
from pygame.locals import *
import numpy as np
import random, sys
import copy
from constant import *
from solve_potential import *
import matplotlib.backends.backend_agg as agg



class Simulation(object):
	def __init__(self):
		pygame.init()
		self.screen = pygame.display.set_mode(SCREENSIZE)
		self.grid = Grid(10**5)
		self.run()



	def run(self):
		running = True
		self.screen.fill(BLACK)
		while running:
			for event in pygame.event.get():
				if event.type == pygame.QUIT:	running = False
				if event.type == pygame.KEYDOWN:
					running = False
			
			self.grid.particle_to_density()
			fig = plt.figure(figsize = (SCREENSIZE[0]/100.0,SCREENSIZE[1]/100.0), dpi = 100, frameon = False)
			#fig, ax = plt.subplots(1)
			fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
			plt.imshow(self.grid.grid.T/10**8, cmap = 'magma', vmin = 0, vmax = 2)
			plt.axis('off')
			plt.axis('tight')
			canvas = agg.FigureCanvasAgg(fig)
			plt.close()
			canvas.draw()
			renderer = canvas.get_renderer()
			raw_data = renderer.tostring_rgb()
			surf = pygame.image.fromstring(raw_data, canvas.get_width_height(), "RGB")
			self.screen.blit(surf, (0,0))
			self.grid.calculate_potential()
			self.grid.calculate_force()
			self.grid.update()
			
			pygame.display.flip()
			
		pygame.quit()
		sys.exit()


if __name__ == "__main__":
	sicheng = Simulation()

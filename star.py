import pygame
from pygame.locals import *
import random
import math

class Star(object):
	def __init__(self, v):
		self.theta = random.random() * math.pi * 2
		self.cos = math.cos(self.theta)
		self.sin = math.sin(self.theta)
		self.z = random.randint(0,10**5)
		self.s = random.randint(0,2500)
		self.o = [540, 337]
		#self.color = (random.randint(0,255),random.randint(0,255),random.randint(0,255))
		self.color = (255, 255, 255)
		self.v = v
		self.x = 0
		self.y = 0

	def trans(self):
		temp = self.s**2 / math.sqrt(self.s**2 + self.z**2)
		self.x = int(self.cos * temp) + self.o[0]
		self.y = int(self.sin * temp) + self.o[1]

	def update(self):
		self.z -= self.v
		if self.z < 0:		self.z = random.randint(10**4, 10**5)
		self.trans()

class Stars(object):
	def __init__(self, n):
		self.n = n
		self.stars = [Star(500) for i in range(n)]

	def update(self):
		for each in self.stars:
			each.update() 


class Game(object):
	def __init__(self):
		pygame.init()
		self.screen = pygame.display.set_mode((1080, 675), pygame.FULLSCREEN)

		self.stars = Stars(10000)
		self.run()
	
	def star(self):
		pass

	def run(self):
		running = True
		

		while running:
			for event in pygame.event.get():
				if event.type == pygame.QUIT:	running = False
				if event.type == pygame.KEYDOWN:
					if event.key == pygame.K_ESCAPE:	running = False	
			self.screen.fill((0, 0, 0))
			self.stars.update()
			for each in self.stars.stars:
				self.screen.fill(each.color, ([each.x, each.y], (1,1)))


			pygame.display.flip()



if __name__ == "__main__":
	sicheng = Game()

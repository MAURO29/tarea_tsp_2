import math

import numpy as np 

from sympy import *


class tablaDH():

	"""
	En esta clase se declarara la tabla de Denavit-Hartenberg.

	Con los parámetros de theta, d, a y alpha se construirá
	la matriz de transformación de 4x4 para cada sistema de referencia.
	"""

	theta=symbols('theta')
	d=symbols('d')
	a=symbols('a')
	alfa=symbols('alpha')
		
	def __init__(self):
			self.theta1,self.theta2,self.theta3,self.theta4=symbols('th1 th2 th3 th4')
			self.l1, self.l2, self.l3, self.l4=symbols('l1 l2 l3 l4')

	def tabla(self):

		"""
        Tabla de DH expresada en forma simbólica
        """
	   
		self.tabla=np.array([
			[self.theta1, self.l3, self.l1, 0],
			[self.theta2, self.l4, self.l2, 0],
			[self.theta3, 0, 0, 0],
			[0, self.theta4, 0, 0]
			])


	

		





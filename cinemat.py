import math

import numpy as np

from sympy import *

from Tabla_DH import *


class Cinematica(tablaDH):
	'''
	Los parámetros Denavit-Hartenberg  resuelven en robótica la cinemática directa.
    
    Permiten establecer la ubicación de los sistemas de referencia de los eslabones en 
    los sistemas robóticos articulados, con cadenas cinemáticas abiertas.
    Definen las transformaciones relativas entre eslabones con tan solo "cuatro 
    parámetros", siendo éste el número mínimo de parámetros para configuraciones
    genéricas.

	Los 4 parámetros son 
		- Angulo θi: :Angulo de la articulación del eje xi −1 al eje i x respecto 
		del eje zi-1. Acepta datos del tipo
		int, float o double.
		- Distancia en Z: Es la distanica desde el origen del sistema de coordenadas
		(i-1)-ésimo hasta la intersección del eje zi-1 con el xi a lo largo del eje
		zi-1. Para su declaración se usan valores del tipo int, float o double.
		- Distancia en X: Distancia de separación desde la intersección del eje 
		zi-1 con el eje xi hasta el orígen del sistema i-ésimo a lo largo del eje xi 
		(o la distancia más corta entre los ejes zi-1 y zi cuando los ejes de
        articulación son paralelos).
        - Angulo αi: Ángulo de separación del eje zi-1 al eje zi respecto del eje xi.
        Al declararla se deben usar valores del tipo int, float o double.

    Para lograr expresar un punto cualquiera del sistema de coordenadas i-ésimo en 
    términos del sistema de coordenadas anterior se proponen las siguientes 
    transformaciones sucesivas:
	1. Girar respecto del eje zi-1 un ángulo θi para alinear el eje xi −1 con el eje 
	xi.
	2. Trasladar a lo largo del eje zi-1 una distancia de di para llevar en coincidencia
	los ejes xi −1 y xi.
	3. Trasladar a lo largo del eje xi una distancia ai para traer en coincidencia 
	tambien los dos orígenes de los ejes x.
	4. Girar respecto del eje xi un ángulo αi para traer en coincidencia a los sistemas
	de coordenadas.

	La multiplicación matricial de transformaciones resulta en otra matriz 
	de 4x4 llamada "Matriz de transformación homogénea".

	'''
 	


	def compute_dh(self):

		'''
		Método que permite calcular la matriz de transformación homogénea simbolicamente.

		Aplicando las 4 transformaciones de traslación y rotación, determina una 
		matriz de 4x4 cuyos parámetros son:
		theta-angúlo de rotación en radianes del tipo int, float o double.
		distancia en z: Traslación en metros, in, ft, etc; del tipo int, float, double.
		distancia en x: Traslación en metros, in, ft, etc; del tipo int, float, double.
		alpha: Ángulo de rotación en radianes del tipo int, float o double.
		'''
		matriz_theta= np.array([
			[cos(tablaDH.theta),- sin(tablaDH.theta),0,0],
			[sin(tablaDH.theta), cos(tablaDH.theta),0,0],
			[0,0,1,0],
			[0,0,0,1]
			])
		matriz_distanciaZ=np.array([
			[1, 0, 0, 0],
			[0, 1, 0, 0],
			[0, 0, 1, tablaDH.d],
			[0, 0, 0, 1]])
		matriz_distanciaX=np.array([
			[1, 0, 0, tablaDH.a],
			[0, 1, 0, 0],
			[0, 0, 1, 0],
			[0, 0, 0, 1]
			])
		matriz_alpha=np.array([
			[1,0,0,0],
			[0, cos(tablaDH.alfa),-sin(tablaDH.alfa),0],
			[0, sin(tablaDH.alfa), cos(tablaDH.alfa),0],
			[0,0,0,1]
			])
		matrizDH=((matriz_theta.dot(matriz_distanciaZ)).dot(matriz_distanciaX)).dot(matriz_alpha)
		
		print(matrizDH)

	def matrizN(self, grados):
		"""
        Matriz de DH para n grados de libertad
        """
		tabla1=tablaDH()
		tabla1.tabla()
		matrizParcial=np.eye(4)
		for i in range(0,grados):
			costh=cos(tabla1.tabla[i,0])
			sinth=sin(tabla1.tabla[i,0])
			d=tabla1.tabla[i,1]
			a=tabla1.tabla[i,2]
			cosalp=cos(tabla1.tabla[i,3])
			sinalp=sin(tabla1.tabla[i,3])

			MatrizDH=np.array([
				[costh,-sinth*cosalp,sinalp*sinth,a*costh],
				[sinth,cosalp*costh,-sinalp*costh,a*sinth],
				[0,sinalp,cosalp,d],
				[0,0,0,1]
				])
			matrizParcial=matrizParcial.dot(MatrizDH)
			j=i+1
			print("Matriz ",i,"A",j)
			print(MatrizDH)
		for i in range(0,4):
			for j in range(0,4):
				matrizParcial[i,j]=trigsimp(matrizParcial[i,j])

		self.A0n=matrizParcial
		print("La matriz de transformación de la base al efector final 0A4 es:")
		print(self.A0n)












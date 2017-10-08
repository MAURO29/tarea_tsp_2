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
    


    def __init__(self,N):
        self.Mi=[0 for x in range(0,N)]
        self.Ai=[0 for x in range(0,N)]

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
        A0n=np.eye(4)
        tabla1=tablaDH()
        tabla1.tabla()
        
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
            
            self.Mi[i]=MatrizDH
            A0n=A0n.dot(MatrizDH)
            self.Ai[i]=A0n
            j=i+1
            print("Matriz de ",i,"a ",j)
            print(self.Mi[i])
            print("\n")
        for i in range(0,4) :
            for j in range(0,4):
                A0n[i,j]=trigsimp(A0n[i,j])
        
        print("La matriz de transformación de la base al efector final es:")
        print(A0n)
        print("\n\n")

    def compute_jacobian(self, N):
        Zi=[0 for i in range(N)]
        ri=[0 for i in range(N)]
        for x in range(0,N):
            Zi[x]=self.Mi[x][0:3,2:3]
            ri[x]=self.Ai[x][0:3,3:4]
        C1=np.cross(Zi[0],ri[N-1],axis=0)
        a=np.asarray(ri[N-1]-ri[0])
        for i in range(0,3) :
            for j in range(0,1):
                a[i,j]=trigsimp(a[i,j])
        C2=np.cross(Zi[1],a,axis=0)
        a=np.asarray(ri[N-1]-ri[1])
        for i in range(0,3) :
            for j in range(0,1):
                a[i,j]=trigsimp(a[i,j])
        C3=np.cross(Zi[2],a,axis=0)
        C4=np.asarray(Zi[3])
        J=np.array([
            [C1[0,0],C2[0,0],C3[0,0],C4[0,0]],
            [C1[1,0],C2[1,0],C3[1,0],C4[1,0]],
            [C1[2,0],C2[2,0],C3[2,0],C4[2,0]],
            [0,0,0,0],
            [0,0,0,0],
            [1,1,1,0]
        ])
        print("El jacobiano del robot es:")
        print(J)
       











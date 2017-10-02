'''Creado miercoles 30 de agosto del 2017
   Autor: González Vázquez Fabian Mauro'''

print("---------------Temas Selectos de Programación-----------")
print("-------------------------Tarea 1------------------------")
print("\n\n")
print("Ejercicio 1:")

import numpy as np
MatrizA01=np.matrix([[0.45,1.3,0.94,1.23],[0,0.83,0.2,2.37],[0.2,0.65,0.4,0.3],[0,0,0,1]])
MatrizA12=np.matrix([[1,0.2,0.85,2.467],[0.54,1.3,0.,0.77],[0.12,0.68,1,0.8],[0,0,0,1]])
vectorx2=np.matrix([[2.3],[0],[24.4],[1]])
print("Matriz A01:")
print(MatrizA01)
print("Matriz A12")
print(MatrizA12)
print("vector x2")
print(vectorx2)
Vectorx0=MatrizA01*MatrizA12*vectorx2
print("el vector resultante de la multiplicación de matrices es:")
print("x0=",Vectorx0)


from cinemat import * 
import numpy as np

'''
Creación del objeto
datos de inciso a)
'''
Matriz_transformacion_homogenea= Cinematica()
'''ángulos de rotación theta'''
theta1=0
theta2=90
theta3=-90
'''traslaciones'''
l1=1
l2=1
l3=2
l4=1
theta4=2

T01=Matriz_transformacion_homogenea.compute_dh(theta1,l3+l4,l1,0)
T12=Matriz_transformacion_homogenea.compute_dh(theta2,0,l2,0)
T23=Matriz_transformacion_homogenea.compute_dh(theta3,theta4,0,0)

vector_origen=np.array([[0],[0],[0],[1]])

print("El producto de matrices de transformación es:")
print(T01@T12@T23)
print("El vector de posición del efector final es:")
print(T01@T12@T23@vector_origen)
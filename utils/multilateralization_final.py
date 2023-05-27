from this import d
import numpy as np
from math import sqrt
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

distances = []
refPoints = []
coordsUnknownPoint = []
# distances = [3.857, 3.988, 3.497] # Replica paper
# coordsUnknownPoint = [26.7726, -1.3389] # Replica paper
# refPoints = [[27.297, -4.953, 1.47], [20.693, -4.849, 1.93], [22.59, 0.524, 1.2], [17.113, -3.003, 2.17], [22.554, 4.727, 1.77]] # Replica paper
N1 = []
N2 = []
t1 = any
t2 = any
A_segundoGrado = 0
B_segundoGrado = 0
C_segundoGrado = 0
x_max = 0
x_min = 0
y_max = 0
y_min = 0

# Ejemplos
# refPoints = [[1, 5.2, 0.3, 9],[5, -3.6, -3, -1],[9, -10, 1, 3], [-2, 3.4, 0.1, 6.1], [-2.9, 10, -5.1, 3.9], [8, 5, -4, 0.7], [-1.4, -2, 3, 6.2]]
# refPoints = [[3.12, -2.2, 1.3, -0.9],[-2.5, 2.6, 6.3, 1.9],[-5, 7.1, -2.1, 5.3], [1.1, 4.3, -5.2, 1.06]]
refPoints = [[1, 5.2, 0.3, 9],[5, -3.6, -3, -1],[9, -10, 1, 3], [-2, 3.4, 0.1, 6.1]]
# refPoints = [[1, 5.2, 0.3, 9],[5, -3.6, -3, -1],[9, -10, 1, 3]]
coordsUnknownPoint = [7, 10, 4.1, 5.7]

# Matriz Identidad
I = np.identity(len(refPoints[0]))

# Matrices para el sist. de ecuaciones
A = np.array([])
b = np.array([])

# Función que calcule las distancias
def calculateDistances():
  global distances
  for coords in refPoints:
    value = 0
    for i in range(len(coordsUnknownPoint)):
      value += (coords[i]-coordsUnknownPoint[i])**2
    s = sqrt(value)
    # s = sqrt((coords[i]-coordsUnknownPoint[i])**2 + (coords[1]-coordsUnknownPoint[1])**2)
    distances.append(s)
      
  print(distances)

# Generar las ecuaciones de A en función de los puntos de referencia
def generateAEcs():
  global A
  for coords in refPoints:
    newRow = [1]
    for value in coords:
      newRow.append(-2*value)
    print(newRow)
    A = np.append(A, np.array(newRow))
    A = np.reshape(A, (-1, len(newRow)))
  print("A =", A)

# Generar las ecuaciones de B en función de los puntos de referencia y distancias
def generateBEcs():
  global b
  for coords, dist in zip(refPoints, distances):
    newRow = dist**2
    for value in coords:
      newRow = newRow - value**2
    b = np.append(b, newRow)
  print("b =", b)

calculateDistances()
generateAEcs()
generateBEcs()

# Calculate EC 2º Grado values
def secondGradeVariables(x_p, x_h):
  global A_segundoGrado, B_segundoGrado, C_segundoGrado

  if len(x_p) != len(x_h):
    raise ValueError("La longitud debe ser la misma")
  else:
    length = len(x_p)

  A_segundoGrado = sum([x_h[i]**2 for i in range(1,length)])
  print(A_segundoGrado)

  B_segundoGrado = sum([2*x_p[i]*x_h[i] for i in range(1,length)]) - x_h[0]
  print(B_segundoGrado)

  C_segundoGrado = sum([x_p[i]**2 for i in range(1,length)]) - x_p[0]
  print(C_segundoGrado)


# Ec 2º Grado
def solve_quadratic(a, b):
  global t1, t2

  solucion1 = (-b-cmath.sqrt(insideSqrt))/(2*a)
  solucion2 = (-b+cmath.sqrt(insideSqrt))/(2*a)

  if insideSqrt < 0:
    t1 = solucion1.real
    t2 = solucion2.real
  else:
    t1 = solucion1.real
    t2 = solucion2.real

def calculateDifferences(x1, x2):
  global d1, d2

  d1 = x1[0] - sum([x1[i]**2 for i in range(1,len(x1))])
  d2 = x2[0] - sum([x2[i]**2 for i in range(1,len(x2))])

  print("dif1 =", d1)
  print("dif2 =", d2)

def findSolution():
  X = np.array(N1)
  Y = np.array(N2)
  P = np.array(refPoints)

  res1 = 0
  res2 = 0
  res = []

  for i in range(len(P)):
    res1 += (np.linalg.norm(X - P[i])**2 - distances[i])**2
    res2 += (np.linalg.norm(Y - P[i])**2 - distances[i])**2

  res.append(res1)
  res.append(res2)

  indexMinValue = np.argmin(res)

  print("minimice =", res[indexMinValue])

  if indexMinValue == 0:
    print("La solución es N1")
  else:
    print("La solución es N2")


# Fin de la definición de funciones, aquí empieza el cálculo del punto desconocido

# x_p (soluciones particulares) y x_h (soluciones homogéneas) usando la pseudo-inversa
A_pseudo_inv = np.linalg.pinv(A)
print(A_pseudo_inv)
x_p = A_pseudo_inv.dot(b)

# SVD (Descomposición en valores singulares). U matriz unitaria. S matriz diagonal. VT matriz transpuesta conjugada.
U, S, VT = np.linalg.svd(A)
x_h = VT[-1] # x_h último vector en la matriz VT (matriz de vectores propios transpuestos)

print("x_p =", x_p)
print("x_h =", x_h)

# Ecuación de 2º Grado
secondGradeVariables(x_p, x_h)
insideSqrt = (B_segundoGrado**2)-(4*A_segundoGrado*C_segundoGrado)
print("insideSqrt =", insideSqrt)

solve_quadratic(A_segundoGrado, B_segundoGrado)

print("t1 =", t1)
print("t2 =", t2)

sol1 = x_p + t1*x_h
sol2 = x_p + t2*x_h

if len(sol1) != len(refPoints[0])+1:
  raise ValueError("Se lía")

calculateDifferences(sol1, sol2)
sol1_copy = sol1
sol2_copy = sol2

sol1 = sol1[1:] # Me interesa quedarme con las coordenadas desconocidas, ya que las primeras lo son.
sol2 = sol2[1:]

print(sol1)
print(sol2)

N1 = sol1.dot(I)
N2 = sol2.dot(I)
print("N1 =", N1)
print("N2 =", N2)

findSolution()

# Representación gráfica

# Calcular los límites del gráfico

def generatePlots(refPoints):
  global x_max, x_min, y_max, y_min, y_max

  resultado = []
  x_values = []
  y_values = []

  for subarray in refPoints:
      resultado.append(subarray[:2])

  for i, (x, y) in enumerate(resultado):
    # if i > 2:
    #   break
    x_values.append(resultado[i][0] + distances[i])
    x_values.append(resultado[i][0] - distances[i])
    y_values.append(resultado[i][1] + distances[i])
    y_values.append(resultado[i][1] - distances[i])
    circle = plt.Circle((x, y), distances[i], color='b', fill=False)
    ax.add_patch(circle)
  
  x_max = max(x_values)
  x_min = min(x_values)
  y_max = max(y_values)
  y_min = min(y_values)
  
# Crear el gráfico
generatePlots(refPoints)
fig, ax = plt.subplots()

# Añadir las circunferencias al gráfico

plt.plot(N1[0], N1[1], 'ro')
plt.plot(N2[0], N2[1], 'yo')
plt.xlabel('x')
plt.ylabel('y')

# Configurar los ejes
print("x_min es:", x_min, "x_max es:", x_max)
print("y_min es:", y_min, "y_max es:", y_max)
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
plt.axis('equal')
plt.show()
## PARA LEER LA PRIMERA FILA DEL CSV QUE ME HA PASADO MATEO.
## LA IDEA ES SIMULAR EL FLUJO DE LA APP, DE LA QUE SABRÍAMOS 2 DE SUS COORDENADAS
## EL OBJETIVO ES CALCULAR LAS COORDENADAS QUE FALTAN PARA PROBAR EL MÉTODO DE MULTILATERALIZACIÓN

import numpy as np
import pandas as pd
from google.colab import files

# Cargar el csv
uploaded = files.upload()

# Lee el csv y extrae los datos
df = pd.read_csv(next(iter(uploaded)))
datos = df.iloc[0].tolist()

# Extrae el nombre del archivo y las coordenadas
nombreFichero = datos[0]
coordenadas_reales_csv = np.array(datos[1:], dtype=float)

# Valores antes de normalizarlos
print("Nombre del fichero:", nombreFichero)
print("Coordenadas:", coordenadas_reales_csv)

# Normaliza los valores
coordenadas_reales_normalizadas_csv = coordenadas_reales_csv / np.max(coordenadas_reales_csv)

# Valores normalizados
print("\nValores normalizados:")
print("Valores normalizados:", coordenadas_reales_normalizadas_csv)

# Array con los dos primeros valores normalizados
coordenadas_conocidas_puntoDesconocido = coordenadas_reales_normalizadas_csv[:2]

# Dos primeros valores normalizados
print("\nDos primeros valores normalizados:")
print("Dos primeros valores normalizados:", coordenadas_conocidas_puntoDesconocido)

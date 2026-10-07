# -*- coding: utf-8 -*-
"""
combinar_carpetas.py
====================
Une los galfitm_resultados.csv de varias carpetas de trabajo (galaxias_1_9,
galaxias_10_22, ...) en una sola tabla, con dos columnas nuevas:

  carpeta       nombre de la carpeta de origen (p. ej. 'galaxias_10_22')
  ruta_carpeta  ruta completa; los scripts la usan para encontrar los FITS
                de cada galaxia (cortes de imagen, apertura, asimetria, bulbo+disco)

Despues, todo el analisis se corre UNA vez sobre la muestra completa, desde la
carpeta donde esta este script (p. ej. analisis_conjunto/):

    python3 combinar_carpetas.py
    python3 procesar_muestra.py
    python3 corregir_ext_kcor.py
    python3 clasificar_morfologia.py
    python3 seleccionar_discos.py
    python3 generar_bd_feedmes.py

Requisito: en cada carpeta de origen ya se corrio GalfitM y leer_output_galfitm.py
(que crea galfitm_resultados.csv).

Uso:
    python3 combinar_carpetas.py                    # busca ../galaxias_*
    python3 combinar_carpetas.py ../galaxias_1_9 ../galaxias_10_22
"""

import os
import sys
import glob
import numpy as np
import pandas as pd

BASE   = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, 'galfitm_resultados.csv')

carpetas = sys.argv[1:] or sorted(set(
    glob.glob(os.path.join(BASE, '..', 'galaxias_*')) + glob.glob(os.path.join(BASE, 'galaxias_*'))))
carpetas = sorted(set(os.path.abspath(c) for c in carpetas if os.path.isdir(c)))
if not carpetas:
    sys.exit('No encontre carpetas galaxias_*. Pasalas como argumento.')

tablas = []
print(f'{"carpeta":25s} {"galaxias":>9s} {"ok":>7s}')
for c in carpetas:
    csv = os.path.join(c, 'galfitm_resultados.csv')
    nombre = os.path.basename(c)
    if not os.path.exists(csv):
        print(f'{nombre:25s}   FALTA galfitm_resultados.csv (corre leer_output_galfitm.py ahi)')
        continue
    t = pd.read_csv(csv, dtype={'objectId': str}, low_memory=False)
    t.insert(1, 'carpeta', nombre)
    t.insert(2, 'ruta_carpeta', c)
    tablas.append(t)
    print(f'{nombre:25s} {len(t):9d} {(t.estado == "ok").sum():7d}')

if not tablas:
    sys.exit('Ninguna carpeta tiene galfitm_resultados.csv')
df = pd.concat(tablas, ignore_index=True)

# Una galaxia en dos carpetas: se conserva la que tiene ajuste ok y menor chi2
dup = df.objectId.duplicated(keep=False)
if dup.any():
    n_ids = df.loc[dup, 'objectId'].nunique()
    df['_orden'] = np.where(df.estado == 'ok', 0, 1)
    df = (df.sort_values(['objectId', '_orden', 'chi2nu'])
            .drop_duplicates('objectId', keep='first')
            .drop(columns='_orden')
            .sort_values(['carpeta', 'objectId'])
            .reset_index(drop=True))
    print(f'\nAVISO: {n_ids} objectId estaban en mas de una carpeta; '
          f'se conservo el ajuste ok con menor chi2')

df.to_csv(SALIDA, index=False)
print(f'\nTabla combinada: {SALIDA}')
print(f'  {len(df)} galaxias de {len(tablas)} carpetas '
      f'({(df.estado == "ok").sum()} con ajuste ok)')
print('\nSiguiente paso:  python3 procesar_muestra.py')

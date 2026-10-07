# -*- coding: utf-8 -*-
"""
filtrar_gri_local.py
====================
Aplica a los lotes ya descargados (v7) el mismo filtro que el notebook v8:
cada galaxia debe tener las bandas g, r e i, y ninguna de ellas puede tener
mas del 30% de pixeles vacios (0 o NaN) en su imagen sci.

Uso (desde la carpeta que contiene galfitm_inputs/):
    python3 filtrar_gri_local.py            # solo lista -> sin_gri.txt
    python3 filtrar_gri_local.py --borrar   # borra las carpetas de la lista
"""

import os
import sys
import shutil
import numpy as np
from astropy.io import fits

DIR = 'galfitm_inputs'
OBLIGATORIAS = ['g', 'r', 'i']
FRAC_VACIA_MAX = 0.30
BORRAR = len(sys.argv) > 1 and sys.argv[1] == '--borrar'

if not os.path.isdir(DIR):
    sys.exit(f'No encuentro la carpeta {DIR}/ (corre el script desde galaxias_1_9)')

malas, revisadas = [], 0
for oid in sorted(os.listdir(DIR)):
    d = os.path.join(DIR, oid)
    if not (os.path.isdir(d) and oid.isdigit()):
        continue
    revisadas += 1
    motivo = []
    for b in OBLIGATORIAS:
        f = os.path.join(d, f'{oid}_sci_{b}.fits')
        if not os.path.exists(f):
            motivo.append(f'sin_{b}')
            continue
        try:
            img = fits.getdata(f)
            if np.mean(~np.isfinite(img) | (img == 0)) > FRAC_VACIA_MAX:
                motivo.append(f'vacia_{b}')
        except Exception:
            motivo.append(f'ilegible_{b}')
    if motivo:
        malas.append((oid, ','.join(motivo)))

with open('sin_gri.txt', 'w') as f:
    for oid, m in malas:
        f.write(f'{oid}\t{m}\n')

print(f'Galaxias revisadas          : {revisadas}')
print(f'Sin g, r, i utiles          : {len(malas)}  (lista en sin_gri.txt)')
print(f'Se conservan                : {revisadas - len(malas)}')

if BORRAR:
    for oid, _ in malas:
        shutil.rmtree(os.path.join(DIR, oid), ignore_errors=True)
    print(f'Borradas {len(malas)} carpetas.')
else:
    print('\nRevisa sin_gri.txt y, si esta bien, corre:')
    print('    python3 filtrar_gri_local.py --borrar')

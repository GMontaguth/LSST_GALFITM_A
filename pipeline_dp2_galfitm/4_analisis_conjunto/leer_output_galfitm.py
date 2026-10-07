# -*- coding: utf-8 -*-
"""
leer_output_galfitm.py
======================
Lee los outputs de GalfitM de TODAS las galaxias en galfitm_inputs/
y genera una tabla con los parametros Sersic por banda.

Estructura esperada (todo dentro de galaxias_1_9):
    galaxias_1_9/
        leer_output_galfitm.py        <- este script
        galfitm_inputs/
            <ID>/
                <ID>_galfitm_out.fits
                <ID>_galfit.01.band   (opcional, para chi2)
                ...

Genera (dentro de galaxias_1_9):
    galfitm_resultados.csv   -> una fila por galaxia (incluye las que fallaron)
    galfitm_imagenes/        -> PNG input/modelo/residuo (si MAKE_IMAGES=True)

Uso:
    python3 leer_output_galfitm.py
"""

import os
import glob
import random
import numpy as np
import pandas as pd
from astropy.io import fits

# =============================================================================
# CONFIG
# =============================================================================
BASE        = os.path.dirname(os.path.abspath(__file__))
GALFITM_DIR = os.path.join(BASE, 'galfitm_inputs')
OUT_CSV     = os.path.join(BASE, 'galfitm_resultados.csv')
OUT_IMG_DIR = os.path.join(BASE, 'galfitm_imagenes')

MAKE_IMAGES = True     # False = mucho mas rapido (solo tabla)
N_IMAGENES  = 100      # cuantas PNG hacer (elegidas al azar). None = todas
SEMILLA     = 42       # misma semilla = mismas galaxias elegidas cada vez

# Catalogo opcional para agregar ra, dec, etc. (deja None si no tienes)
CATALOGO    = None     # ej: os.path.join(BASE, 'dp2_galaxias.csv')
COLS_CAT    = ['ra', 'dec']   # columnas del catalogo a copiar

BANDAS     = ['u', 'g', 'r', 'i', 'z', 'y']
LAMBDA_EFF = {'u': 3550, 'g': 4670, 'r': 6160,
              'i': 7470, 'z': 8920, 'y': 10200}

# Parametros Sersic (componente 1) y cielo (componente 2)
PARAMS_SERSIC = ['XC', 'YC', 'MAG', 'RE', 'N', 'AR', 'PA']

if MAKE_IMAGES:
    import matplotlib
    matplotlib.use('Agg')           # no abre ventanas
    import matplotlib.pyplot as plt
    os.makedirs(OUT_IMG_DIR, exist_ok=True)


# =============================================================================
# HELPERS
# =============================================================================

def parsear_valor(raw):
    """
    Parsea valores del header de GalfitM:
      '300.8425 +/- 0.0078'  -> (300.8425, 0.0078, flag=0)
      '[301.0000]'           -> (301.0, 0.0, flag=0)   parametro fijo
      '*19.5 +/- 0.3*'       -> (19.5, 0.3, flag=1)    GalfitM marca problema
    """
    raw    = str(raw).strip()
    flag   = 1 if '*' in raw else 0
    fijo   = raw.startswith('[') and raw.endswith(']')
    limpio = raw.replace('*', '').strip('[]() ')

    if '+/-' in limpio:
        a, b = limpio.split('+/-', 1)
        try:
            return float(a), (0.0 if fijo else float(b)), flag
        except ValueError:
            return np.nan, np.nan, flag
    try:
        return float(limpio), (0.0 if fijo else np.nan), flag
    except ValueError:
        return np.nan, np.nan, flag


def header_parametros(hdul):
    """Busca el HDU cuyo header contiene los parametros (1_MAG_*)."""
    for hdu in hdul:
        if any(k.startswith('1_MAG_') for k in hdu.header.keys()):
            return hdu.header
    return None


def get_chi2nu(hdul, carpeta, oid):
    for hdu in hdul:
        if 'CHI2NU' in hdu.header:
            try:
                return float(hdu.header['CHI2NU'])
            except (ValueError, TypeError):
                pass
    # Respaldo: archivo .band (el script de corrida lo renombra como ID_galfit.01.band)
    for band_path in glob.glob(os.path.join(carpeta, '*galfit.*.band')):
        try:
            with open(band_path) as f:
                for line in f:
                    if 'Chi^2/nu' in line:
                        return float(line.split('=')[1].split(',')[0].strip())
        except (OSError, ValueError, IndexError):
            pass
    return np.nan


def leer_parametros(header, banda):
    fila = {}
    bu = banda.upper()
    flags = 0
    for p in PARAMS_SERSIC:
        val, err, flag = np.nan, np.nan, 0
        key = f'1_{p}_{bu}'
        if key in header:
            val, err, flag = parsear_valor(header[key])
        fila[f'{p}_{banda}']   = val
        fila[f'e_{p}_{banda}'] = err
        flags += flag
    # Cielo (componente 2)
    key = f'2_SKY_{bu}'
    fila[f'SKY_{banda}'] = parsear_valor(header[key])[0] if key in header else np.nan
    fila[f'flag_{banda}'] = flags   # >0: GalfitM marco algun parametro con *
    return fila


def asinh_norm(img, lo=0.5, hi=99.5, stretch=0.5):
    v1, v2 = np.nanpercentile(img, [lo, hi])
    x = np.clip((img - v1) / (v2 - v1 + 1e-12), 0, 1)
    return np.clip(np.arcsinh(x / stretch) / np.arcsinh(1.0 / stretch), 0, 1)


def hacer_imagen(hdul, oid, chi2nu, bandas):
    ext = {hdu.name.upper(): hdu for hdu in hdul if hdu.data is not None}
    bandas = [b for b in bandas if all(f'{t}_{b.upper()}' in ext
                                       for t in ('INPUT', 'MODEL', 'RESIDUAL'))]
    if not bandas:
        return None

    nb = len(bandas)
    fig, axes = plt.subplots(3, nb, figsize=(nb * 3.2, 9.5),
                             facecolor='black', squeeze=False)
    for col, b in enumerate(bandas):
        for row, tipo in enumerate(['INPUT', 'MODEL', 'RESIDUAL']):
            ax   = axes[row, col]
            data = ext[f'{tipo}_{b.upper()}'].data.astype(np.float32)
            if tipo == 'RESIDUAL':
                vmax = np.nanpercentile(np.abs(data), 99) or 1.0
                ax.imshow(data, cmap='RdBu_r', origin='lower', vmin=-vmax, vmax=vmax)
            else:
                ax.imshow(asinh_norm(data), cmap='gray', origin='lower', vmin=0, vmax=1)
            if row == 0:
                ax.set_title(f'{b}  {LAMBDA_EFF.get(b, 0)} A', color='white', fontsize=9)
            if col == 0:
                ax.set_ylabel(tipo, color='white', fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_facecolor('black')

    chi = f'{chi2nu:.2f}' if np.isfinite(chi2nu) else 'N/A'
    fig.suptitle(f'objectId={oid}  |  Chi2/nu={chi}', color='white',
                 fontsize=10, fontweight='bold')
    plt.tight_layout()
    out_png = os.path.join(OUT_IMG_DIR, f'{oid}_galfitm.png')
    plt.savefig(out_png, dpi=110, facecolor='black')
    plt.close(fig)
    return out_png


# =============================================================================
# MAIN
# =============================================================================

cat_dict = {}
if CATALOGO and os.path.exists(CATALOGO):
    cat = pd.read_parquet(CATALOGO) if CATALOGO.endswith('.parquet') else pd.read_csv(CATALOGO)
    cat['objectId'] = cat['objectId'].astype(str)
    cat_dict = cat.set_index('objectId').to_dict('index')
    print(f'Catalogo cargado: {len(cat)} objetos')
else:
    COLS_CAT = []   # sin catalogo no se agregan columnas vacias

carpetas = sorted(d for d in os.listdir(GALFITM_DIR)
                  if os.path.isdir(os.path.join(GALFITM_DIR, d)))
total = len(carpetas)
print(f'Galaxias en {GALFITM_DIR}: {total}')

# Elegir al azar las galaxias (con output) que tendran imagen
con_output = [d for d in carpetas
              if os.path.exists(os.path.join(GALFITM_DIR, d, f'{d}_galfitm_out.fits'))]
if MAKE_IMAGES and N_IMAGENES is not None and N_IMAGENES < len(con_output):
    random.seed(SEMILLA)
    elegidas_img = set(random.sample(con_output, N_IMAGENES))
else:
    elegidas_img = set(con_output)
if MAKE_IMAGES:
    print(f'Se haran {len(elegidas_img)} imagenes PNG (al azar)')
print()

resultados = []
conteo = {'ok': 0, 'sin_output': 0, 'error': 0}
mostrado_debug = False

for i, oid in enumerate(carpetas, 1):
    carpeta   = os.path.join(GALFITM_DIR, oid)
    fits_path = os.path.join(carpeta, f'{oid}_galfitm_out.fits')

    fila = {'objectId': oid, 'estado': '', 'chi2nu': np.nan,
            'n_bandas': 0, 'bandas': ''}
    for c in COLS_CAT:
        fila[c] = cat_dict.get(oid, {}).get(c, np.nan)

    if not os.path.exists(fits_path) or os.path.getsize(fits_path) < 1000:
        fila['estado'] = 'sin_output'
        conteo['sin_output'] += 1
        resultados.append(fila)
        continue

    try:
        with fits.open(fits_path) as hdul:
            if not mostrado_debug:
                print('=== Estructura del primer output ===')
                for j, hdu in enumerate(hdul):
                    print(f'  HDU {j}: {hdu.name}')
                mostrado_debug = True

            header = header_parametros(hdul)
            if header is None:
                raise ValueError('no se encontraron parametros (1_MAG_*) en el header')

            nombres = [hdu.name.upper() for hdu in hdul]
            bandas  = [b for b in BANDAS if f'INPUT_{b.upper()}' in nombres
                       or f'1_MAG_{b.upper()}' in header]

            fila['chi2nu']   = get_chi2nu(hdul, carpeta, oid)
            fila['n_bandas'] = len(bandas)
            fila['bandas']   = ','.join(bandas)
            for b in bandas:
                fila.update(leer_parametros(header, b))
            fila['estado'] = 'ok'

            if MAKE_IMAGES and oid in elegidas_img:
                hacer_imagen(hdul, oid, fila['chi2nu'], bandas)

        conteo['ok'] += 1
        if i % 100 == 0 or i == total:
            print(f'  [{i}/{total}] procesadas...')

    except Exception as e:
        fila['estado'] = f'error: {e}'
        conteo['error'] += 1
        print(f'  ERROR {oid}: {e}')

    resultados.append(fila)


# =============================================================================
# GUARDAR
# =============================================================================
df = pd.DataFrame(resultados)

cols_info = ['objectId', 'estado', 'chi2nu', 'n_bandas', 'bandas'] + COLS_CAT
# Ordena parametros por banda (i, z, y...) y luego por parametro
cols_par = []
for b in BANDAS:
    for p in PARAMS_SERSIC:
        cols_par += [f'{p}_{b}', f'e_{p}_{b}']
    cols_par += [f'SKY_{b}', f'flag_{b}']
cols_par = [c for c in cols_par if c in df.columns]
df = df[cols_info + cols_par]

df.to_csv(OUT_CSV, index=False)

print()
print('=' * 60)
print(f'  OK          : {conteo["ok"]}')
print(f'  Sin output  : {conteo["sin_output"]}')
print(f'  Con error   : {conteo["error"]}')
print(f'  Tabla       : {OUT_CSV}')
if MAKE_IMAGES:
    print(f'  Imagenes    : {OUT_IMG_DIR}/')
print('=' * 60)

prev = ['objectId', 'chi2nu', 'MAG_i', 'e_MAG_i', 'RE_i', 'N_i', 'AR_i', 'PA_i', 'flag_i']
print(df[df.estado == 'ok'][[c for c in prev if c in df.columns]].head(10).to_string(index=False))

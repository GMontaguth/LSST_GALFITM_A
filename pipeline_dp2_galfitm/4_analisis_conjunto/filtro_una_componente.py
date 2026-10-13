# -*- coding: utf-8 -*-
"""
filtro_una_componente.py
========================
Segundo paso del filtro bulbo + disco: detecta automáticamente las galaxias que pasan
filtro_bd.py (bd_ok, bulbo_no_resuelto) pero que en realidad se describen con UNA sola
componente, y las reasigna a la categoría que ya usa el resto del pipeline.

Criterios (banda de referencia BANDA, parámetros de GalfitM):

  1. SIN CRUCE -> 'esferoide_puro'
     El bulbo domina el brillo superficial en TODOS los radios entre 0.5 FWHM y
     R_MAX_RED x R_e del disco: el "disco" nunca domina, solo corrige el perfil del
     bulbo (típico de elípticas ajustadas con B+D). Perfiles 1D a lo largo del eje
     mayor de cada componente, calculados con su magnitud, R_e, n, b/a:
        Sérsic: I(r) = I_e exp{-b_n [(r/R_e)^(1/n) - 1]}
        exponencial: I(r) = I_0 exp(-r/R_s)
     Estas galaxias quedan fuera del análisis de bulbos y discos (entran a las
     Partes 1 y 2 con su ajuste de una componente).

  2. DISCO DUPLICADO -> 'disco_puro'
     Las dos componentes son el mismo disco repartido en dos exponenciales:
        n_b < DUP_N_MAX,  R_e,b / R_e,d > DUP_RATIO_MIN,
        |(b/a)_b - (b/a)_d| < DUP_DQ_MAX,  |PA_b - PA_d| < DUP_DPA_MAX
        (el PA no se exige si las dos son casi redondas, b/a > DUP_Q_REDONDA)
     Pasan a 'disco_puro': masa_tamano_componentes.py ya usa para ellas la magnitud y
     el R_e del ajuste de una componente como disco, y no les asigna bulbo.

Columnas nuevas en galfitm_bd_filtrado.csv:
  categoria_bd_antes   categoría que había puesto filtro_bd.py
  una_comp_motivo      'sin_cruce', 'disco_duplicado' o vacío
  corte_sin_cruce, corte_disco_duplicado   (booleanas, para todas las galaxias evaluadas)
  r_disco_domina_arcsec  radio desde el que domina el disco (NaN si nunca)

También reescribe muestra_bd_ok.csv y deja revisar_sin_cruce.csv y
revisar_disco_duplicado.csv (20 al azar) para verificarlas a ojo con
  python3 inspeccion_visual_bd.py --ids analisis_galfitm/bd/filtro/revisar_sin_cruce.csv --incluir-auto

Orden del pipeline (desde analisis_conjunto/):
  python3 leer_output_bd.py
  python3 filtro_bd.py
  python3 filtro_una_componente.py          <- este
  python3 masa_tamano_componentes.py
  (inspección visual y notebook)
Si vuelves a correr filtro_bd.py, vuelve a correr este después.

Uso:
  python3 filtro_una_componente.py
  python3 filtro_una_componente.py --solo-ver      # no modifica nada, solo cuenta
"""
import os
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy.special import gammaincinv, gamma
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIG
# =============================================================================
BASE = os.path.dirname(os.path.abspath(__file__))
AG = os.path.join(BASE, 'analisis_galfitm')
FILTRO = os.path.join(AG, 'bd', 'filtro', 'galfitm_bd_filtrado.csv')
BD_OK = os.path.join(AG, 'bd', 'filtro', 'muestra_bd_ok.csv')
RESULT = os.path.join(AG, 'bd', 'galfitm_bd_resultados.csv')
MORFO = os.path.join(AG, 'morfologia', 'galaxias_morfologia.csv')
OUT_DIR = os.path.join(AG, 'bd', 'filtro')

BANDA = 'r'                    # si falta, se usa BANDA_RESPALDO
BANDA_RESPALDO = 'i'
PIXSCALE = 0.2
APLICAR_A = ['bd_ok', 'bulbo_no_resuelto']
FWHM_DEFECTO_PIX = 3.5

# 1. sin cruce
R_MAX_RED = 3.0                # se mira hasta 3 R_e del disco

# 2. disco duplicado
DUP_N_MAX = 1.5
DUP_RATIO_MIN = 0.6
DUP_DQ_MAX = 0.1
DUP_DPA_MAX = 10.0
DUP_Q_REDONDA = 0.85

N_REVISAR = 20
SEMILLA = 42

ap = argparse.ArgumentParser()
ap.add_argument('--solo-ver', action='store_true', help='no escribe nada, solo cuenta')
args = ap.parse_args()


def leer(p):
    if not os.path.exists(p):
        sys.exit(f'No encuentro {p}')
    return pd.read_csv(p, dtype={'objectId': str}, low_memory=False)


f = leer(FILTRO)
if 'categoria_bd_antes' in f:            # ya se corrió antes: partir de lo que dejó filtro_bd.py
    f['categoria_bd'] = f['categoria_bd_antes']
else:
    f['categoria_bd_antes'] = f['categoria_bd']


# -----------------------------------------------------------------------------
# parámetros: primero los de la tabla filtrada (son los del ajuste elegido);
# los que falten se toman de galfitm_bd_resultados.csv
# -----------------------------------------------------------------------------
def nombres(b):
    return {'B_MAG': f'B_MAG_{b}', 'B_RE': f'B_RE_{b}', 'B_N': f'B_N_{b}', 'B_AR': f'B_AR_{b}',
            'B_PA': f'B_PA_{b}', 'D_MAG': f'D_MAG_{b}', 'D_RS': f'D_RS_{b}', 'D_RE': f'D_RE_{b}',
            'D_AR': f'D_AR_{b}', 'D_PA': f'D_PA_{b}'}


faltan = [c for b in (BANDA, BANDA_RESPALDO) for c in nombres(b).values() if c not in f]
if faltan and os.path.exists(RESULT):
    r = leer(RESULT)
    llave = ['objectId', 'ajuste'] if ('ajuste' in f and 'ajuste' in r) else ['objectId']
    traer = [c for c in faltan if c in r]
    if traer:
        r = r[llave + traer].drop_duplicates(llave)
        f = f.merge(r, on=llave, how='left')


def col(c):
    return pd.to_numeric(f[c], errors='coerce') if c in f else pd.Series(np.nan, index=f.index)


P = {}
for k in nombres(BANDA):
    P[k] = col(nombres(BANDA)[k]).fillna(col(nombres(BANDA_RESPALDO)[k]))
P['D_RS'] = P['D_RS'].fillna(P['D_RE'] / 1.678)
P['D_RE'] = P['D_RE'].fillna(1.678 * P['D_RS'])
# FWHM de la PSF en píxeles: FWHM_PSF primero. RE_B_FWHM (cociente R_e,b / FWHM) y FWHM_KPC no son la FWHM.
c_fwhm = [c for c in f.columns if c.lower().startswith('fwhm') and 'kpc' not in c.lower()]
FWHM = pd.Series(np.nan, index=f.index)
for c in sorted(c_fwhm, key=lambda c: (c != 'FWHM_PSF', c)):
    FWHM = FWHM.fillna(pd.to_numeric(f[c], errors='coerce'))
FWHM = FWHM.fillna(FWHM_DEFECTO_PIX)
print(f'FWHM de la PSF: {"columna " + ", ".join(c_fwhm) if c_fwhm else f"valor por defecto {FWHM_DEFECTO_PIX} px"}')

evaluar = f.categoria_bd.isin(APLICAR_A).to_numpy()
completo = np.all([np.isfinite(P[k]) for k in ('B_MAG', 'B_RE', 'B_N', 'D_MAG', 'D_RS')], axis=0)
print(f'Galaxias en {APLICAR_A}: {int(evaluar.sum())}  (con parámetros completos: {int((evaluar & completo).sum())})')


# -----------------------------------------------------------------------------
# 1. sin cruce
# -----------------------------------------------------------------------------
def perfiles(i, r):
    Fb = 10 ** (-0.4 * P['B_MAG'].iat[i])
    Fd = 10 ** (-0.4 * P['D_MAG'].iat[i])
    n = np.clip(P['B_N'].iat[i], 0.2, 10)
    re, rs = P['B_RE'].iat[i], P['D_RS'].iat[i]
    qb = P['B_AR'].iat[i] if np.isfinite(P['B_AR'].iat[i]) else 1.0
    qd = P['D_AR'].iat[i] if np.isfinite(P['D_AR'].iat[i]) else 1.0
    bn = gammaincinv(2 * n, 0.5)
    Ie = Fb / (2 * np.pi * qb * re ** 2 * np.exp(bn) * n * bn ** (-2 * n) * gamma(2 * n))
    I0 = Fd / (2 * np.pi * qd * rs ** 2)
    Ib = Ie * np.exp(-bn * ((r / re) ** (1 / n) - 1))
    Id = I0 * np.exp(-r / rs)
    return Ib, Id


sin_cruce = np.zeros(len(f), bool)
r_domina = np.full(len(f), np.nan)
for i in np.where(evaluar & completo)[0]:
    r_lo = 0.5 * FWHM.iat[i]
    r_hi = max(R_MAX_RED * P['D_RE'].iat[i], 2 * r_lo)
    r = np.geomspace(r_lo, r_hi, 300)
    Ib, Id = perfiles(i, r)
    domina = Id > Ib
    if domina.any():
        r_domina[i] = r[np.argmax(domina)] * PIXSCALE
    else:
        sin_cruce[i] = True

# -----------------------------------------------------------------------------
# 2. disco duplicado
# -----------------------------------------------------------------------------
dpa = np.abs(((P['B_PA'] - P['D_PA']) + 90) % 180 - 90)
redondas = (P['B_AR'] > DUP_Q_REDONDA) & (P['D_AR'] > DUP_Q_REDONDA)
duplicado = (evaluar & completo
             & (P['B_N'] < DUP_N_MAX)
             & (P['B_RE'] / P['D_RE'] > DUP_RATIO_MIN)
             & ((P['B_AR'] - P['D_AR']).abs() < DUP_DQ_MAX)
             & ((dpa < DUP_DPA_MAX) | redondas)).to_numpy().copy()
duplicado &= ~sin_cruce            # si cumple las dos, manda la de esferoide

f['corte_sin_cruce'] = sin_cruce
f['corte_disco_duplicado'] = duplicado
f['r_disco_domina_arcsec'] = r_domina
f['una_comp_motivo'] = np.select([sin_cruce, duplicado], ['sin_cruce', 'disco_duplicado'], '')
f.loc[sin_cruce, 'categoria_bd'] = 'esferoide_puro'
f.loc[duplicado, 'categoria_bd'] = 'disco_puro'

# -----------------------------------------------------------------------------
# resumen
# -----------------------------------------------------------------------------
print('\nReasignadas:')
print(f'  sin cruce (el disco nunca domina)  -> esferoide_puro : {int(sin_cruce.sum())}')
print(f'  disco duplicado                    -> disco_puro     : {int(duplicado.sum())}')
print('\nCategoría antes -> después (solo las que cambiaron):')
cambio = f.categoria_bd != f.categoria_bd_antes
if cambio.any():
    print(pd.crosstab(f.categoria_bd_antes[cambio], f.una_comp_motivo[cambio]).to_string())
print('\nCategorías finales:')
print(f.categoria_bd.value_counts().to_string())
tipo = f['TypeName'] if 'TypeName' in f else None
if tipo is None and os.path.exists(MORFO):
    mo = leer(MORFO)
    if 'TypeName' in mo:
        tipo = f[['objectId']].merge(mo[['objectId', 'TypeName']].drop_duplicates('objectId'),
                                      on='objectId', how='left').TypeName.set_axis(f.index)
if tipo is not None and cambio.any():
    print('\nReasignadas por tipo morfológico:')
    print(pd.crosstab(f.una_comp_motivo[cambio], tipo[cambio]).to_string())

if args.solo_ver:
    print('\n(--solo-ver: no se escribió nada)')
    sys.exit(0)

f.to_csv(FILTRO, index=False)
f[f.categoria_bd == 'bd_ok'].to_csv(BD_OK, index=False)
rng = np.random.default_rng(SEMILLA)
for motivo in ('sin_cruce', 'disco_duplicado'):
    q = f[f.una_comp_motivo == motivo]
    if len(q):
        q.sample(min(N_REVISAR, len(q)), random_state=SEMILLA)[['objectId', 'categoria_bd', 'una_comp_motivo']] \
            .to_csv(os.path.join(OUT_DIR, f'revisar_{motivo}.csv'), index=False)
print(f'\nActualizado: {FILTRO}')
print(f'             {BD_OK}  ({int((f.categoria_bd == "bd_ok").sum())} bd_ok)')
print('Para verificar a ojo una muestra:')
print('  python3 inspeccion_visual_bd.py --ids analisis_galfitm/bd/filtro/revisar_sin_cruce.csv --incluir-auto')
print('  python3 inspeccion_visual_bd.py --ids analisis_galfitm/bd/filtro/revisar_disco_duplicado.csv --incluir-auto')
print('Después corre masa_tamano_componentes.py para que las masas y tamaños usen las categorías nuevas.')

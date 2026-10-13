# -*- coding: utf-8 -*-
"""
filtro_bd.py
============
Filtro logico de las descomposiciones bulbo + disco (GalfitM, receta Vika 2014)
pensado para la pregunta: ¿el tamano del bulbo sigue la relacion masa-tamano?
(Lima-Dias et al. 2024 vs Mendez-Abreu et al. 2021).

Cada galaxia recibe UNA categoria (se aplica en este orden; la primera que
cumple gana) y ademas columnas booleanas corte_* con cada criterio por separado:

  fallida            ajuste no confiable: Re, n, b/a del bulbo o Rs, b/a del disco
                     marcados con '*' por GalfitM, chi2/nu fuera de [CHI2_MIN, CHI2_MAX],
                     |m_tot(B+D) - m_single-Sersic| > DMAG_TOT_MAX, o B/T sin valor
  perturbada         asimetria de statmorph A > ASIM_MAX (0.35, Conselice 2003): irregular,
                     interaccion, vecina o estrella encima -> el modelo B+D no aplica
  invertida          Re_bulbo >= RATIO_RE_MAX * Re_disco (el "bulbo" es la componente
                     extendida: descomposicion invertida, candidata a re-ajuste)
  esferoide_puro     B/T >= BT_MAX  (el disco no aporta: no hay dos componentes)
  disco_puro         B/T <= BT_MIN  (el bulbo no aporta)
  disco_absurdo      Re_disco > DISCO_RE_MAX * Re_single-Sersic (disco que modela cielo/halo)
  posible_barra      n_b < BARRA_N_MAX, (b/a)_b < BARRA_Q_MAX, |PA_b - PA_d| > BARRA_DPA_MIN,
                     Re_b < Re_d  (componente central plana y desalineada: barra, no bulbo;
                     se excluye para comparar con Mendez-Abreu 2021, muestra sin barras)
  bulbo_aguja        (b/a)_bulbo < B_AR_MIN (0.2): barra, disco interno o artefacto
  disco_de_canto     (b/a)_disco < D_AR_MIN (0.3): galaxia de canto, B+D no confiable
  bulbo_absurdo      n_b < N_B_MIN o n_b > N_B_MAX
  bulbo_no_resuelto  Re_b < FWHM_MIN * FWHM_PSF  (el tamano del bulbo no se puede medir:
                     estos bulbos aplanan artificialmente la relacion masa-tamano)
  bd_ok              pasa todo  -> muestra para el analisis

La FWHM de la PSF se mide en cada galaxia (<ID>_psf_<b>.fits) y se guarda en
psf_fwhm.csv (solo se miden las nuevas).

Entradas (analisis_galfitm/ junto al script):
  bd/galfitm_bd_resultados.csv            (leer_output_bd.py)
  muestra/muestra_clasificada.csv         (procesar_muestra.py; single-Sersic + categoria)
     o crossmatch/galfitm_x_lsst.csv
  morfologia/galaxias_morfologia.csv      (opcional: TypeName, logM_taylor, z)

Salidas (analisis_galfitm/bd/filtro/):
  galfitm_bd_filtrado.csv      todas, con categoria_bd y columnas corte_*
  muestra_bd_ok.csv            solo bd_ok
  revisar_<categoria>.csv      20 al azar de cada categoria (para mirar las figuras:
                               python3 leer_output_bd.py --ids .../revisar_<cat>.csv)
  psf_fwhm.csv                 FWHM de la PSF por galaxia y banda
  muestra_componentes.csv      bd_ok + esferoides puros + discos puros, con columnas
                               comunes: RE_ESF, N_ESF, RE_DISCO, MAG_ESF_<b>, MAG_DISCO_<b>
                               (en las puras = parametros single-Sersic) y grupo_comp
  22_filtro_bd.png             diagnosticos
  23_componentes_masa_re.png   Re de esferoides y discos vs M* de la galaxia:
                               B+D confiables vs puras

Uso (desde analisis_conjunto/):
    python3 filtro_bd.py
    python3 filtro_bd.py --ba-min 0.5        # ademas exige (b/a)_disco >= 0.5 en bd_ok
    python3 filtro_bd.py --fwhm-min 1.5      # resolucion mas exigente
    (las variantes se guardan en filtro/ba_min_0.5/, filtro/fwhm_min_1.5/, ...;
     psf_fwhm.csv se comparte)
"""
import os
import sys
import glob
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIG
# =============================================================================
BASE = os.path.dirname(os.path.abspath(__file__))
AG = os.path.join(BASE, 'analisis_galfitm')
BD_CSV = os.path.join(AG, 'bd', 'galfitm_bd_resultados.csv')
BD2_CSV = os.path.join(AG, 'bd', 'galfitm_bd2_resultados.csv')   # 2do ajuste (si existe)
# 2do ajuste: valores pegados al borde de sus restricciones (generar_bd_refit.py) = no confiable
LIMITE_BD2 = dict(ratio_max=0.88, n_min=0.32, n_max=7.8, dmag_max=2.95)
RANGO_CAT = {'bd_ok': 0, 'disco_inclinado': 1, 'bulbo_no_resuelto': 2, 'disco_puro': 3,
             'esferoide_puro': 4, 'posible_barra': 5, 'bulbo_aguja': 5, 'disco_de_canto': 5,
             'disco_absurdo': 6, 'perturbada': 6,
             'bulbo_absurdo': 7, 'invertida': 8, 'fallida': 9}
SS_CSVS = [os.path.join(AG, 'muestra', 'muestra_clasificada.csv'),
           os.path.join(AG, 'crossmatch', 'galfitm_x_lsst.csv')]
MORFO_CSV = os.path.join(AG, 'morfologia', 'galaxias_morfologia.csv')
OUT_DIR = os.path.join(AG, 'bd', 'filtro')
PSF_CACHE = os.path.join(OUT_DIR, 'psf_fwhm.csv')

BANDAS_COMP = ['g', 'r', 'i']   # magnitudes por componente (para masas: g-i y M_i)
BANDA = 'r'                 # banda de referencia (si falta en una galaxia se usa BANDA_ALT)
BANDA_ALT = 'i'
PIXSCALE = 0.2              # arcsec/px
FWHM_DEFECTO = 3.5          # px (~0.7"), si no se encuentra la PSF

CHI2_MIN, CHI2_MAX = 0.5, 5.0
DMAG_TOT_MAX = 0.3          # |m(B+D) - m(single Sersic)| en la banda de referencia
RATIO_RE_MAX = 1.0          # Re_b / Re_d
BT_MIN, BT_MAX = 0.06, 0.94
DISCO_RE_MAX = 3.0          # Re_d / Re_ss
N_B_MIN, N_B_MAX = 0.3, 8.0
FWHM_MIN = 1.0              # Re_b >= FWHM_MIN * FWHM_PSF
BARRA_N_MAX, BARRA_Q_MAX, BARRA_DPA_MIN = 0.7, 0.5, 20.0

# parametros estructurales: si GalfitM marca alguno con '*' el ajuste no sirve
PARAMS_CRITICOS = ('B_RE', 'B_N', 'B_AR', 'D_RS', 'D_RE', 'D_AR')

# Controles agregados tras revisar las figuras del 2do ajuste:
B_AR_MIN = 0.2       # bulbo con b/a < 0.2 = "aguja": barra, disco interno o artefacto
D_AR_MIN = 0.3       # disco con b/a < 0.3 = galaxia de canto: B+D no confiable
# Asimetria de statmorph (Conselice 2003: A > 0.35 = fusion / interaccion). Detecta
# irregulares, galaxias en interaccion y mezclas con vecinas o estrellas brillantes.
STATMORPH_CSV = os.path.join(AG, 'statmorph_resultados.csv')
ASIM_MAX = 0.35
ASIM_BANDAS = ('r', 'i', 'g')    # primera banda disponible con flag <= 1 y S/N >= 2.5
SM_FLAG_MAX, SM_SN_MIN = 1, 2.5

N_REVISAR = 20
SEMILLA = 42

CATEGORIAS = ['fallida', 'perturbada', 'invertida', 'esferoide_puro', 'disco_puro',
              'disco_absurdo', 'posible_barra', 'bulbo_aguja', 'disco_de_canto',
              'bulbo_absurdo', 'bulbo_no_resuelto', 'bd_ok']
NOMBRE_EN = {'fallida': 'failed fit', 'invertida': 'inverted (Re_b ≥ Re_d)',
             'esferoide_puro': f'pure spheroid (B/T ≥ {BT_MAX})',
             'disco_puro': f'pure disc (B/T ≤ {BT_MIN})', 'disco_absurdo': 'oversized disc',
             'posible_barra': 'possible bar', 'bulbo_absurdo': 'unphysical bulge n',
             'bulbo_no_resuelto': 'unresolved bulge', 'bd_ok': 'reliable B+D',
             'perturbada': f'disturbed (A > {ASIM_MAX})',
             'bulbo_aguja': f'needle-like bulge (b/a < {B_AR_MIN})',
             'disco_de_canto': f'edge-on disc (b/a < {D_AR_MIN})'}
COLOR_CAT = {'fallida': '#7f8c8d', 'invertida': '#8e44ad', 'esferoide_puro': '#c0392b',
             'disco_puro': '#2e86c1', 'disco_absurdo': '#16a085', 'posible_barra': '#f39c12',
             'bulbo_absurdo': '#d35400', 'bulbo_no_resuelto': '#bdc3c7', 'bd_ok': '#27ae60',
             'perturbada': '#34495e', 'bulbo_aguja': '#e84393', 'disco_de_canto': '#00b894'}

ap = argparse.ArgumentParser()
ap.add_argument('--fwhm-min', type=float, default=FWHM_MIN,
                help=f'Re_b minimo en unidades de FWHM de la PSF (defecto {FWHM_MIN})')
ap.add_argument('--ba-min', type=float, default=None,
                help='exigir ademas (b/a) del disco >= este valor en bd_ok')
ap.add_argument('--banda', default=BANDA, help=f'banda de referencia (defecto {BANDA})')
args = ap.parse_args()
FWHM_MIN, BANDA = args.fwhm_min, args.banda
os.makedirs(OUT_DIR, exist_ok=True)
if args.ba_min is not None or args.fwhm_min != 1.0:      # variantes: subcarpeta propia
    sufijo = []
    if args.ba_min is not None:
        sufijo.append(f'ba_min_{args.ba_min}')
    if args.fwhm_min != 1.0:
        sufijo.append(f'fwhm_min_{args.fwhm_min}')
    OUT_DIR = os.path.join(OUT_DIR, '_'.join(sufijo))
    os.makedirs(OUT_DIR, exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.3})


# =============================================================================
# 1. LEER Y UNIR
# =============================================================================
if not os.path.exists(BD_CSV):
    sys.exit(f'No existe {BD_CSV}\nCorre antes: python3 leer_output_bd.py')
bd = pd.read_csv(BD_CSV, dtype={'objectId': str}, low_memory=False).assign(ajuste='bd')
print(f'Bulbo+disco leidas: {len(bd)}')
if os.path.exists(BD2_CSV):
    bd2 = pd.read_csv(BD2_CSV, dtype={'objectId': str}, low_memory=False).assign(ajuste='bd2')
    print(f'  2do ajuste (_bd2): {len(bd2)}  -> se elige el mejor de los dos por galaxia')
    bd = pd.concat([bd, bd2], ignore_index=True)

ss_csv = next((p for p in SS_CSVS if os.path.exists(p)), None)
if ss_csv is None:
    sys.exit('No encuentro la tabla single-Sersic (muestra_clasificada.csv o galfitm_x_lsst.csv)')
ss = pd.read_csv(ss_csv, dtype={'objectId': str}, low_memory=False)
cols_ss = ['objectId'] + [c for c in ss.columns
                          if c.split('_')[0] in ('MAG', 'RE', 'N', 'AR', 'PA')
                          and c.split('_')[-1] in set(BANDAS_COMP) | {BANDA, BANDA_ALT}
                          and c.count('_') == 1]
if 'categoria' in ss:
    cols_ss.append('categoria')
ss = ss[cols_ss].drop_duplicates('objectId').rename(
    columns={c: f'SS_{c}' for c in cols_ss if c not in ('objectId', 'categoria')})
df = bd.merge(ss, on='objectId', how='left')
print(f'  con single-Sersic ({os.path.basename(ss_csv)}): {df[f"SS_MAG_{BANDA}"].notna().sum()}')

if os.path.exists(MORFO_CSV):
    mo = pd.read_csv(MORFO_CSV, dtype={'objectId': str}, low_memory=False)
    keep = ['objectId'] + [c for c in ('TypeName', 'logM_taylor', 'z') if c in mo]
    df = df.merge(mo[keep].drop_duplicates('objectId'), on='objectId', how='left')
    print(f'  con morfologia/masa: {df.get("logM_taylor", pd.Series(dtype=float)).notna().sum()}')


# =============================================================================
# 2. BANDA DE REFERENCIA POR GALAXIA
# =============================================================================
# Re, n, b/a y PA son constantes en lambda (receta Vika), asi que basta una banda;
# se usa BANDA y, si la galaxia no la tiene, BANDA_ALT.
def col(p, b):
    return df[p + '_' + b] if p + '_' + b in df else pd.Series(np.nan, index=df.index)

usa_alt = col('BT', BANDA).isna() & col('BT', BANDA_ALT).notna()
df['banda_ref'] = np.where(usa_alt, BANDA_ALT, BANDA)
for p in ('B_MAG', 'B_RE', 'B_N', 'B_AR', 'B_PA', 'D_MAG', 'D_RE', 'D_AR', 'D_PA',
          'BT', 'MAG_TOT', 'RATIO_RE', 'SS_MAG', 'SS_RE', 'SS_N', 'SS_AR'):
    df[p] = col(p, BANDA).where(~usa_alt, col(p, BANDA_ALT))


# =============================================================================
# 3. FWHM DE LA PSF (medida en cada galaxia, con cache)
# =============================================================================
def fwhm_psf(path):
    """FWHM = 2 sqrt(A/pi), A = area de pixeles >= max/2 (robusto a PSF no gaussianas)."""
    try:
        p = fits.getdata(path).astype(float)
        p = np.nan_to_num(p)
        return 2.0 * np.sqrt(np.sum(p >= 0.5 * p.max()) / np.pi)
    except Exception:
        return np.nan


raices = sorted({os.path.abspath(c) for c in
                 glob.glob(os.path.join(BASE, 'galaxias_*')) +
                 glob.glob(os.path.join(BASE, '..', 'galaxias_*'))})
ruta_carpeta = {os.path.basename(c): c for c in raices}

cache = (pd.read_csv(PSF_CACHE, dtype={'objectId': str})
         if os.path.exists(PSF_CACHE) else pd.DataFrame(columns=['objectId', 'banda', 'fwhm_px']))
cache = cache[cache.fwhm_px.notna()]                    # las que no tenian PSF se reintentan
hechas = set(zip(cache.objectId, cache.banda))
nuevas = []
pend = df[[(o, b) not in hechas for o, b in zip(df.objectId, df.banda_ref)]]
if len(pend):
    print(f'Midiendo FWHM de la PSF en {len(pend)} galaxias...')
for k, (o, b, c) in enumerate(zip(pend.objectId, pend.banda_ref, pend.carpeta), 1):
    raiz = ruta_carpeta.get(str(c))
    f = os.path.join(raiz, 'galfitm_inputs', o, f'{o}_psf_{b}.fits') if raiz else ''
    nuevas.append({'objectId': o, 'banda': b,
                   'fwhm_px': fwhm_psf(f) if f and os.path.exists(f) else np.nan})
    if k % 500 == 0:
        print(f'  [{k}/{len(pend)}]')
if nuevas:
    cache = pd.concat([cache, pd.DataFrame(nuevas)], ignore_index=True)
    cache.to_csv(PSF_CACHE, index=False)
fw = cache.drop_duplicates(['objectId', 'banda'], keep='last')
df = df.merge(fw.rename(columns={'banda': 'banda_ref', 'fwhm_px': 'FWHM_PSF'}),
              on=['objectId', 'banda_ref'], how='left')
sin_psf = df.FWHM_PSF.isna().sum()
med_fwhm = df.FWHM_PSF.median()
df['FWHM_PSF'] = df.FWHM_PSF.fillna(med_fwhm if np.isfinite(med_fwhm) else FWHM_DEFECTO)
print(f'FWHM PSF mediana: {df.FWHM_PSF.median():.2f} px = {df.FWHM_PSF.median()*PIXSCALE:.2f}"'
      + (f'  ({sin_psf} sin PSF: se usa la mediana)' if sin_psf else ''))


# =============================================================================
# 4. CORTES
# =============================================================================
flags = df.get('params_flag', pd.Series('', index=df.index)).fillna('').astype(str)
df['corte_flag_estructural'] = flags.apply(
    lambda s: any(x.rsplit('_', 1)[0] in PARAMS_CRITICOS for x in s.split(';') if x))
df['corte_chi2'] = ~df.chi2nu.between(CHI2_MIN, CHI2_MAX)
df['DMAG_TOT'] = df.MAG_TOT - df.SS_MAG
df['corte_dmag'] = df.DMAG_TOT.abs() > DMAG_TOT_MAX            # NaN (sin SS) no corta
df['corte_bt_nan'] = df[['BT', 'B_RE', 'B_N', 'D_RE']].isna().any(axis=1)   # sin parametros
df['corte_invertida'] = df.RATIO_RE >= RATIO_RE_MAX
df['corte_esferoide'] = df.BT >= BT_MAX
df['corte_disco_puro'] = df.BT <= BT_MIN
df['corte_disco_absurdo'] = df.D_RE > DISCO_RE_MAX * df.SS_RE
dpa = (df.B_PA - df.D_PA).abs() % 180
df['DPA'] = np.minimum(dpa, 180 - dpa)
df['corte_barra'] = ((df.B_N < BARRA_N_MAX) & (df.B_AR < BARRA_Q_MAX) &
                     (df.DPA > BARRA_DPA_MIN) & (df.RATIO_RE < 1))
df['corte_n_bulbo'] = (df.B_N < N_B_MIN) | (df.B_N > N_B_MAX)
df['RE_B_FWHM'] = df.B_RE / df.FWHM_PSF
df['corte_resolucion'] = df.RE_B_FWHM < FWHM_MIN

# --- asimetria de statmorph
df['ASIM'] = np.nan
if os.path.exists(STATMORPH_CSV):
    sm = pd.read_csv(STATMORPH_CSV, dtype={'objectId': str}, low_memory=False)
    if {'band', 'asymmetry'} <= set(sm.columns):
        ok_sm = pd.to_numeric(sm.get('flag'), errors='coerce') <= SM_FLAG_MAX
        if 'sn_per_pixel' in sm:
            ok_sm &= pd.to_numeric(sm.sn_per_pixel, errors='coerce') >= SM_SN_MIN
        sm = sm[ok_sm]
        for b in ASIM_BANDAS[::-1]:                 # la primera de la lista tiene prioridad
            a = sm[sm.band == b].drop_duplicates('objectId', keep='last').set_index('objectId').asymmetry
            v = df.objectId.map(a)
            df['ASIM'] = v.where(v.notna(), df.ASIM)
    print(f'Asimetria de statmorph: {df.loc[df.ASIM.notna(), "objectId"].nunique()} de '
          f'{df.objectId.nunique()} galaxias')
else:
    print(f'AVISO: no existe {STATMORPH_CSV}: sin corte de asimetria')
df['corte_asimetria'] = df.ASIM > ASIM_MAX                     # sin statmorph: no corta
df['corte_bulbo_aguja'] = df.B_AR < B_AR_MIN
df['corte_de_canto'] = df.D_AR < D_AR_MIN

es2 = df.ajuste == 'bd2'
df['corte_en_limite'] = es2 & ((df.RATIO_RE > LIMITE_BD2['ratio_max']) |
                               (df.B_N < LIMITE_BD2['n_min']) | (df.B_N > LIMITE_BD2['n_max']) |
                               ((df.B_MAG - df.D_MAG).abs() > LIMITE_BD2['dmag_max']))
orden = [('fallida', df.corte_flag_estructural | df.corte_chi2 | df.corte_dmag | df.corte_bt_nan
          | df.corte_en_limite),
         ('perturbada', df.corte_asimetria),
         ('invertida', df.corte_invertida), ('esferoide_puro', df.corte_esferoide),
         ('disco_puro', df.corte_disco_puro), ('disco_absurdo', df.corte_disco_absurdo),
         ('posible_barra', df.corte_barra), ('bulbo_aguja', df.corte_bulbo_aguja),
         ('disco_de_canto', df.corte_de_canto), ('bulbo_absurdo', df.corte_n_bulbo),
         ('bulbo_no_resuelto', df.corte_resolucion)]
df['categoria_bd'] = np.select([m.fillna(False).astype(bool) for _, m in orden],
                               [c for c, _ in orden], default='bd_ok')
if args.ba_min is not None:
    df['corte_ba_disco'] = df.D_AR < args.ba_min
    df.loc[(df.categoria_bd == 'bd_ok') & df.corte_ba_disco, 'categoria_bd'] = 'disco_inclinado'
    CATEGORIAS.insert(-1, 'disco_inclinado')
    NOMBRE_EN['disco_inclinado'] = f'inclined disc (b/a < {args.ba_min})'
    COLOR_CAT['disco_inclinado'] = '#95a5a6'

# Con 2do ajuste: una fila por galaxia, la de mejor categoria (empate -> el 1er ajuste)
rescate = None
if df.ajuste.nunique() > 1:
    antes = df[df.ajuste == 'bd'].set_index('objectId').categoria_bd
    df['_r'] = df.categoria_bd.map(RANGO_CAT).fillna(9)
    df['_a'] = (df.ajuste == 'bd2').astype(int)
    df = (df.sort_values(['objectId', '_r', '_a']).drop_duplicates('objectId')
            .drop(columns=['_r', '_a']).reset_index(drop=True))
    usa2 = df[df.ajuste == 'bd2']
    rescate = pd.crosstab(antes.reindex(usa2.objectId).fillna('sin_1er_ajuste').values,
                          usa2.categoria_bd.values, rownames=['1er ajuste'], colnames=['2do ajuste'])

# =============================================================================
# 4b. TABLA DE COMPONENTES: bulbo+disco confiables + puras (single-Sersic)
# =============================================================================
# esferoide = bulbo (bd_ok) o la galaxia entera (esferoide_puro, parametros single-Sersic)
# disco     = disco (bd_ok) o la galaxia entera (disco_puro, parametros single-Sersic)
GRUPO = {'bd_ok': 'bulge+disc', 'esferoide_puro': 'pure spheroid', 'disco_puro': 'pure disc'}
df['grupo_comp'] = df.categoria_bd.map(GRUPO).fillna('')
es_bd, es_esf, es_dis = (df.categoria_bd == k for k in ('bd_ok', 'esferoide_puro', 'disco_puro'))
nan = pd.Series(np.nan, index=df.index)
for p_bd, p_ss, dest in (('B_RE', 'SS_RE', 'RE_ESF'), ('B_N', 'SS_N', 'N_ESF'),
                         ('B_AR', 'SS_AR', 'AR_ESF')):
    df[dest] = df[p_bd].where(es_bd, df.get(p_ss, nan).where(es_esf))
for p_bd, p_ss, dest in (('D_RE', 'SS_RE', 'RE_DISCO'), ('D_AR', 'SS_AR', 'AR_DISCO')):
    df[dest] = df[p_bd].where(es_bd, df.get(p_ss, nan).where(es_dis))
for b in BANDAS_COMP:
    df[f'MAG_ESF_{b}'] = col('B_MAG', b).where(es_bd, col('SS_MAG', b).where(es_esf))
    df[f'MAG_DISCO_{b}'] = col('D_MAG', b).where(es_bd, col('SS_MAG', b).where(es_dis))
# las puras tambien deben estar resueltas y, si se pide, con disco poco inclinado
df['ESF_RESUELTO'] = df.RE_ESF >= FWHM_MIN * df.FWHM_PSF
if args.ba_min is not None:
    df.loc[es_dis & (df.AR_DISCO < args.ba_min), 'grupo_comp'] = ''

# tamanos fisicos (kpc) si hay z: para la relacion masa-tamano
if 'z' in df:
    try:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
        z = pd.to_numeric(df.z, errors='coerce')
        ok = z > 0
        kpc_arcsec = np.full(len(df), np.nan)
        kpc_arcsec[ok.values] = cosmo.kpc_proper_per_arcmin(z[ok].values).value / 60.0
        df['RE_B_KPC'] = df.B_RE * PIXSCALE * kpc_arcsec
        df['RE_D_KPC'] = df.D_RE * PIXSCALE * kpc_arcsec
        df['FWHM_KPC'] = df.FWHM_PSF * PIXSCALE * kpc_arcsec
        df['RE_ESF_KPC'] = df.RE_ESF * PIXSCALE * kpc_arcsec
        df['RE_DISCO_KPC'] = df.RE_DISCO * PIXSCALE * kpc_arcsec
    except ImportError:
        pass


# =============================================================================
# 5. GUARDAR
# =============================================================================
df.to_csv(os.path.join(OUT_DIR, 'galfitm_bd_filtrado.csv'), index=False)
ok = df[df.categoria_bd == 'bd_ok']
ok.to_csv(os.path.join(OUT_DIR, 'muestra_bd_ok.csv'), index=False)
comp = df[df.grupo_comp != '']
cols_comp = (['objectId', 'carpeta', 'grupo_comp', 'categoria_bd', 'BT', 'FWHM_PSF',
              'RE_ESF', 'N_ESF', 'AR_ESF', 'RE_DISCO', 'AR_DISCO', 'ESF_RESUELTO']
             + [f'MAG_{c}_{b}' for c in ('ESF', 'DISCO') for b in BANDAS_COMP]
             + ['RE_ESF_KPC', 'RE_DISCO_KPC', 'FWHM_KPC', 'SS_N_' + BANDA,
                'TypeName', 'logM_taylor', 'z'])
comp[[c for c in cols_comp if c in comp]].to_csv(
    os.path.join(OUT_DIR, 'muestra_componentes.csv'), index=False)
rng = np.random.default_rng(SEMILLA)
cols_rev = ['objectId', 'carpeta', 'ajuste', 'categoria_bd', 'chi2nu', 'BT', 'RATIO_RE', 'B_RE', 'B_N',
            'B_AR', 'D_RE', 'D_AR', 'DPA', 'ASIM', 'RE_B_FWHM', 'DMAG_TOT', 'params_flag']
for c in CATEGORIAS:
    s = df[df.categoria_bd == c]
    if len(s):
        s = s.iloc[rng.choice(len(s), min(N_REVISAR, len(s)), replace=False)]
        s[[x for x in cols_rev if x in s]].to_csv(
            os.path.join(OUT_DIR, f'revisar_{c}.csv'), index=False)


# =============================================================================
# 6. FIGURA
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))
ax = axes[0, 0]
bins = np.linspace(0, 1, 41)
ax.hist(df.BT.dropna(), bins=bins, histtype='step', color='k', lw=1.5, label=f'all ({df.BT.notna().sum()})')
ax.hist(ok.BT.dropna(), bins=bins, color=COLOR_CAT['bd_ok'], alpha=0.6, label=f'reliable ({len(ok)})')
for v in (BT_MIN, BT_MAX):
    ax.axvline(v, color='r', ls='--')
ax.set_xlabel(f'B/T ({BANDA})'); ax.set_ylabel('N'); ax.legend(fontsize=8)
ax.set_title('Bulge-to-total ratio')

ax = axes[0, 1]
lr = np.log10(df.RATIO_RE.where(df.RATIO_RE > 0))
bins = np.linspace(np.nanpercentile(lr, 0.5), np.nanpercentile(lr, 99.5), 50)
ax.hist(lr.dropna(), bins=bins, histtype='step', color='k', lw=1.5, label='all')
ax.hist(np.log10(ok.RATIO_RE.where(ok.RATIO_RE > 0)).dropna(), bins=bins,
        color=COLOR_CAT['bd_ok'], alpha=0.6, label='reliable')
ax.axvline(np.log10(RATIO_RE_MAX), color='r', ls='--')
ax.set_xlabel(r'log $R_{e,b}/R_{e,d}$'); ax.legend(fontsize=8)
ax.set_title('Bulge / disc size ratio')

ax = axes[0, 2]
ln = np.log10(df.B_N.where(df.B_N > 0))
bins = np.linspace(-1.05, 1.2, 46)
ax.hist(ln.dropna(), bins=bins, histtype='step', color='k', lw=1.5, label='all')
ax.hist(np.log10(ok.B_N.where(ok.B_N > 0)).dropna(), bins=bins, color=COLOR_CAT['bd_ok'],
        alpha=0.6, label='reliable')
for v in (N_B_MIN, N_B_MAX, BARRA_N_MAX):
    ax.axvline(np.log10(v), color='r' if v != BARRA_N_MAX else 'orange', ls='--')
ax.set_xlabel(r'log $n_b$'); ax.legend(fontsize=8); ax.set_title('Bulge Sérsic index')

ax = axes[1, 0]
for c in CATEGORIAS:
    s = df[df.categoria_bd == c]
    ax.scatter(s.RE_B_FWHM, s.B_N, s=4, alpha=0.5, color=COLOR_CAT[c], label=NOMBRE_EN[c])
ax.axvline(FWHM_MIN, color='k', ls='--', lw=1)
ax.axvline(1.5, color='k', ls=':', lw=1)
ax.axhline(N_B_MIN, color='r', ls='--', lw=1)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$R_{e,b}$ / FWHM$_{\rm PSF}$'); ax.set_ylabel(r'$n_b$')
ax.set_title('Bulge size vs PSF and Sérsic index'); ax.legend(fontsize=6, markerscale=3)

ax = axes[1, 1]
s = df[df.categoria_bd.isin(['bd_ok', 'bulbo_no_resuelto', 'invertida'])]
sc = ax.scatter(s.D_RE * PIXSCALE, s.B_RE * PIXSCALE, c=s.BT, s=5, cmap='coolwarm',
                vmin=0, vmax=1, alpha=0.7)
lim = np.nanpercentile(np.r_[s.D_RE, s.B_RE] * PIXSCALE, [0.5, 99.5])
ax.plot(lim, lim, 'k--', lw=1, label='1:1')
ax.axhline(FWHM_MIN * med_fwhm * PIXSCALE if np.isfinite(med_fwhm) else FWHM_DEFECTO * PIXSCALE,
           color='gray', ls=':', label=f'{FWHM_MIN}×FWHM')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$R_{e,d}$ [arcsec]'); ax.set_ylabel(r'$R_{e,b}$ [arcsec]')
ax.set_title('Bulge vs disc size (reliable, unresolved, inverted)')
ax.legend(fontsize=8); plt.colorbar(sc, ax=ax, label='B/T')

ax = axes[1, 2]
cnt = df.categoria_bd.value_counts().reindex(CATEGORIAS).fillna(0)
ax.barh([NOMBRE_EN[c] for c in CATEGORIAS], cnt.values, color=[COLOR_CAT[c] for c in CATEGORIAS])
for i, v in enumerate(cnt.values):
    ax.text(v, i, f' {int(v)} ({100*v/len(df):.0f}%)', va='center', fontsize=8)
ax.invert_yaxis(); ax.set_xlabel('N galaxies'); ax.set_title('Classification (first criterion met)')
ax.set_xlim(0, cnt.max() * 1.3)

fig.suptitle(f'Bulge + disc logical filter  (N = {len(df)}, reliable = {len(ok)})', fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, '22_filtro_bd.png'), dpi=120)
plt.close(fig)


# Figura 23: tamano de esferoides y discos vs masa de la galaxia,
# bulbo+disco confiables comparados con las puras (single-Sersic)
if 'RE_ESF_KPC' in df and 'logM_taylor' in df:
    M = pd.to_numeric(df.logM_taylor, errors='coerce')
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
    paneles = [(axes[0], 'RE_ESF_KPC', [('bulge+disc', 'bulges (B+D)', '#27ae60'),
                                        ('pure spheroid', f'pure spheroids (B/T ≥ {BT_MAX})', '#c0392b')],
                'Spheroidal component'),
               (axes[1], 'RE_DISCO_KPC', [('bulge+disc', 'discs (B+D)', '#2e86c1'),
                                          ('pure disc', f'pure discs (B/T ≤ {BT_MIN})', '#e67e22')],
                'Disc component')]
    bins = np.arange(8.0, 12.01, 0.5)
    for ax, cre, grupos, tit in paneles:
        for g, lab, cc in grupos:
            m = (df.grupo_comp == g) & M.notna() & (df[cre] > 0)
            if cre == 'RE_ESF_KPC':
                m &= df.ESF_RESUELTO
            x, y = M[m], np.log10(df.loc[m, cre])
            ax.scatter(x, y, s=5, alpha=0.35, color=cc, label=f'{lab} (N={m.sum()})')
            xc, yc, e = [], [], []
            for a, b in zip(bins[:-1], bins[1:]):
                s_ = y[(x >= a) & (x < b)]
                if len(s_) >= 10:
                    bs = [np.median(np.random.default_rng(i).choice(s_, len(s_))) for i in range(300)]
                    xc.append((a + b) / 2); yc.append(np.median(s_)); e.append(np.std(bs))
            if xc:
                ax.errorbar(xc, yc, yerr=e, color=cc, marker='o', ms=7, mec='k', lw=2, capsize=3)
        if 'FWHM_KPC' in df:
            ax.axhline(np.log10(FWHM_MIN * df.FWHM_KPC.median()), color='gray', ls=':',
                       label=f'{FWHM_MIN}×PSF FWHM (median)')
        ax.set_xlabel(r'log $M_\star$ (galaxy) [$M_\odot$]')
        ax.set_ylabel(r'log $R_e$ [kpc]'); ax.set_title(tit); ax.legend(fontsize=8)
    fig.suptitle('Component sizes vs galaxy stellar mass: B+D decompositions vs single-component '
                 'galaxies (medians in 0.5 dex bins, N ≥ 10)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, '23_componentes_masa_re.png'), dpi=120)
    plt.close(fig)


# =============================================================================
# 7. RESUMEN
# =============================================================================
print('\n' + '=' * 70)
print(f'  Galaxias bulbo+disco : {len(df)}   (banda {BANDA}; {usa_alt.sum()} usan {BANDA_ALT})')
print('-' * 70)
for c in CATEGORIAS:
    n = (df.categoria_bd == c).sum()
    print(f'  {c:18s} {n:6d}  ({100*n/len(df):4.1f}%)')
if rescate is not None:
    print('-' * 70)
    print(f'  Galaxias que usan el 2do ajuste: {int(rescate.values.sum())}  (categoria 1ro -> 2do):')
    print('  ' + rescate.to_string().replace('\n', '\n  '))
print('-' * 70)
print('  Cada corte por separado (una galaxia puede fallar varios):')
for c in [x for x in df.columns if x.startswith('corte_')]:
    print(f'    {c:26s} {int(df[c].fillna(False).astype(bool).sum()):6d}')
print('-' * 70)
if len(ok):
    print(f'  bd_ok: B/T mediana {ok.BT.median():.2f}, n_b mediana {ok.B_N.median():.2f}, '
          f'Re_b/Re_d mediana {ok.RATIO_RE.median():.2f}')
    print(f'  bd_ok con Re_b >= 1.5 FWHM: {(ok.RE_B_FWHM >= 1.5).sum()}  '
          f'(>= 2 FWHM: {(ok.RE_B_FWHM >= 2).sum()})')
    if 'logM_taylor' in ok:
        m = pd.to_numeric(ok.logM_taylor, errors='coerce')
        print(f'  bd_ok con log M* >= 9.5: {(m >= 9.5).sum()}')
    if 'RE_B_KPC' in ok:
        print(f'  Re_b mediano: {ok.RE_B_KPC.median():.2f} kpc   '
              f'(FWHM PSF mediana: {ok.FWHM_KPC.median():.2f} kpc)')
print('-' * 70)
print('  Tabla de componentes (muestra_componentes.csv):')
for g in ('bulge+disc', 'pure spheroid', 'pure disc'):
    s_ = df[df.grupo_comp == g]
    extra = ''
    if g != 'bulge+disc' and f'SS_N_{BANDA}' in s_:
        extra = f'   n single-Sersic mediana: {s_[f"SS_N_{BANDA}"].median():.2f}'
    print(f'    {g:15s} {len(s_):6d}{extra}')
if 'TypeName' in df:
    print('\n  Categoria vs tipo morfologico:')
    print(pd.crosstab(df.categoria_bd, df.TypeName.fillna('sin_tipo'), margins=True).to_string())
print('=' * 70)
print(f'Resultados en {OUT_DIR}/')
print('Para mirar ejemplos de una categoria:')
print(f'  python3 leer_output_bd.py --ids {os.path.relpath(OUT_DIR, BASE)}/revisar_invertida.csv')

# -*- coding: utf-8 -*-
"""
procesar_muestra.py
===================
Hace en un solo script los pasos 2, 3 y 4 del analisis:

  PASO 2  cortes de calidad de GalfitM      (antes: analizar_galfitm.py)
  PASO 3  crossmatch con el catalogo LSST   (antes: crossmatch_lsst.py)
  PASO 4  muestra confiable A/B/C/D         (antes: muestra_confiable.py)
          (incluye el corte de bandas obligatorias g y r)

Necesita, en la misma carpeta que este script:
  galfitm_resultados.csv     (salida de leer_output_galfitm.py, paso 1)
  galfitm_inputs/            (carpetas de las galaxias con los FITS)
  el parquet de LSST         (ruta en LSST_PARQUET, abajo)

Todo se guarda en analisis_galfitm/ (subcarpetas crossmatch/ y muestra/).
Los umbrales de cada corte estan al inicio de cada paso (busca "CONFIG").

Uso:
    python3 procesar_muestra.py              # corre los pasos 2, 3 y 4
    python3 procesar_muestra.py 3 4          # corre solo los pasos indicados
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
warnings.simplefilter('ignore', pd.errors.PerformanceWarning)
warnings.simplefilter('ignore', RuntimeWarning)

# =============================================================================
# CONFIGURACION GENERAL  (lo que normalmente cambias entre muestras)
# =============================================================================
DIR_MUESTRA = os.path.dirname(os.path.abspath(__file__))
LSST_PARQUET          = os.path.join(DIR_MUESTRA, 'z01_lsst_1.parquet')  # <- ajusta
REVISAR_IMAGENES_FITS = True    # paso 2: medir bandas vacias / fondo malo en los FITS
RECALCULAR_METRICAS   = False   # paso 4: True = rehacer apertura y residuos siempre
                                #   (si las guardadas son de otra muestra se rehacen solas)




# Carpeta de cada galaxia. Si la tabla viene de combinar_carpetas.py tiene la
# columna 'ruta_carpeta' (la galaxias_X_Y de origen); si no, se usa esta carpeta.
def dir_galaxia(oid, ruta=None):
    base = ruta if isinstance(ruta, str) and ruta else DIR_MUESTRA
    return os.path.join(base, 'galfitm_inputs', str(oid))


def hay_imagenes(df):
    if 'ruta_carpeta' in df:
        return df['ruta_carpeta'].dropna().map(os.path.isdir).any()
    return os.path.isdir(os.path.join(DIR_MUESTRA, 'galfitm_inputs'))


# #############################################################################
# PASO 2: CORTES DE CALIDAD DE GALFITM
# #############################################################################
def paso2():

    # =============================================================================
    # CONFIG
    # =============================================================================
    BASE        = os.path.dirname(os.path.abspath(__file__))
    CSV_IN      = os.path.join(BASE, 'galfitm_resultados.csv')
    GALFITM_DIR = os.path.join(BASE, 'galfitm_inputs')
    OUT_DIR     = os.path.join(BASE, 'analisis_galfitm')

    REVISAR_IMAGENES = REVISAR_IMAGENES_FITS  # (configurado arriba)
    PIXSCALE         = 0.2     # arcsec/pixel

    # ---- Cortes por GALAXIA ----
    CHI2_MIN = 0.5             # chi2nu ~ 0 => algo raro (sigma mal, banda vacia...)
    CHI2_MAX = 5.0

    # ---- Cortes por BANDA ----
    MAG_MIN, MAG_MAX = 10.0, 25.0
    RE_MIN,  RE_MAX  = 0.5, 150.0      # pix (imagen de 601 px)
    N_MIN,   N_MAX   = 0.21, 7.9       # cerca de los limites de GalfitM = no confiable
    AR_MIN           = 0.1
    COLOR_NSIGMA     = 5.0             # outlier de color respecto a la poblacion
    COLOR_MIN_DESV   = 1.0             # minimo en mag para llamarlo outlier

    # ---- Cortes de imagen (solo si REVISAR_IMAGENES) ----
    FRAC_VACIA_MAX = 0.3   # fraccion maxima de pixeles = 0 o NaN en INPUT
    FONDO_MAX      = 3.0   # variacion de fondo a gran escala / ruido (gradientes)

    BANDAS     = ['u', 'g', 'r', 'i', 'z', 'y']
    LAMBDA_EFF = {'u': 3550, 'g': 4670, 'r': 6160, 'i': 7470, 'z': 8920, 'y': 10200}
    COLORES    = {'u': '#7b2cbf', 'g': '#2a9d8f', 'r': '#e76f51',
                  'i': '#d62828', 'z': '#6d3a1f', 'y': '#333333'}
    PARAMS     = ['XC', 'YC', 'MAG', 'RE', 'N', 'AR', 'PA']

    os.makedirs(OUT_DIR, exist_ok=True)


    # =============================================================================
    # 1. CALIDAD DE IMAGENES (opcional)
    # =============================================================================
    def metricas_imagen(img, bloque=20):
        """
        frac_vacia : fraccion de pixeles 0 o NaN
        fondo      : (p90 - p10 de medianas en bloques) / ruido
                     ruido medido con diferencias entre pixeles vecinos (no le
                     afectan los gradientes). Un fondo normal da ~0.3-1;
                     un gradiente o imagen suave sin ruido da valores grandes.
        """
        img = np.asarray(img, dtype=float)
        vacia = ~np.isfinite(img) | (img == 0)
        frac_vacia = vacia.mean()
        if frac_vacia > 0.95:
            return frac_vacia, np.nan

        x = np.where(vacia, np.nan, img)
        dif = np.diff(x, axis=1).ravel()
        dif = dif[np.isfinite(dif)]
        ruido = 1.4826 * np.median(np.abs(dif - np.median(dif))) / np.sqrt(2)

        ny, nx = (s // bloque * bloque for s in x.shape)
        bl = x[:ny, :nx].reshape(ny // bloque, bloque, nx // bloque, bloque)
        med = np.nanmedian(bl, axis=(1, 3)).ravel()
        med = med[np.isfinite(med)]
        if ruido <= 0 or len(med) < 10:
            return frac_vacia, np.inf
        return frac_vacia, (np.percentile(med, 90) - np.percentile(med, 10)) / ruido


    def revisar_imagenes(ids, rutas=None):
        from astropy.io import fits
        rutas = rutas or {}
        filas = []
        for k, oid in enumerate(ids, 1):
            path = os.path.join(dir_galaxia(oid, rutas.get(oid)), f'{oid}_galfitm_out.fits')
            fila = {'objectId': oid}
            try:
                with fits.open(path) as hdul:
                    nombres = {h.name.upper(): h for h in hdul}
                    for b in BANDAS:
                        h = nombres.get(f'INPUT_{b.upper()}')
                        if h is not None and h.data is not None:
                            fv, fo = metricas_imagen(h.data)
                            fila[f'frac_vacia_{b}'] = fv
                            fila[f'fondo_{b}'] = fo
            except Exception as e:
                print(f'  [img] {oid}: {e}')
            filas.append(fila)
            if k % 500 == 0:
                print(f'  imagenes revisadas: {k}/{len(ids)}')
        return pd.DataFrame(filas)


    # =============================================================================
    # 2. LEER Y APLICAR CORTES
    # =============================================================================
    df = pd.read_csv(CSV_IN, dtype={'objectId': str})
    print(f'Galaxias en tabla: {len(df)}')

    # --- calidad de imagen
    if REVISAR_IMAGENES and hay_imagenes(df):
        print('Revisando imagenes INPUT (puede tardar unos minutos)...')
        ok_ = df[df.estado == 'ok']
        rutas = (dict(zip(ok_.objectId, ok_.ruta_carpeta)) if 'ruta_carpeta' in ok_ else None)
        # Cache incremental: solo se miden las galaxias nuevas (p. ej. al sumar una carpeta)
        cache_q = os.path.join(OUT_DIR, 'calidad_imagenes.csv')
        q = (pd.read_csv(cache_q, dtype={'objectId': str})
             if os.path.exists(cache_q) else pd.DataFrame({'objectId': pd.Series([], dtype=str)}))
        nuevas = [o for o in ok_['objectId'] if o not in set(q.objectId)]
        if nuevas:
            print(f'  {len(q)} ya medidas, {len(nuevas)} nuevas')
            q = pd.concat([q, revisar_imagenes(nuevas, rutas)], ignore_index=True)
            q.to_csv(cache_q, index=False)
        else:
            print(f'  todas las imagenes ya estaban medidas ({cache_q})')
        df = df.merge(q, on='objectId', how='left')
    elif REVISAR_IMAGENES:
        print('No encuentro las carpetas galfitm_inputs: se omite la revision de imagenes')

    # --- cortes por galaxia
    motivo_gal = pd.Series('', index=df.index)
    motivo_gal[df.estado != 'ok'] = 'sin_output'
    motivo_gal[(df.estado == 'ok') & (df.chi2nu < CHI2_MIN)] = 'chi2_bajo'
    motivo_gal[(df.estado == 'ok') & (df.chi2nu > CHI2_MAX)] = 'chi2_alto'
    df['galaxia_ok'] = motivo_gal == ''
    df['motivo_galaxia'] = motivo_gal

    # --- outliers de color: MAG_b - mediana(MAG de la galaxia), comparado con la poblacion
    mags = df[[f'MAG_{b}' for b in BANDAS if f'MAG_{b}' in df]]
    med_gal = mags.where((mags > MAG_MIN) & (mags < MAG_MAX)).median(axis=1)

    df_clean = df.copy()
    motivos_banda = {oid: [] for oid in df.objectId}

    for b in BANDAS:
        if f'MAG_{b}' not in df:
            continue
        tiene = df[f'MAG_{b}'].notna()
        malo = pd.Series('', index=df.index)

        def marcar(cond, txt):
            sel = tiene & cond & (malo == '')
            malo[sel] = txt

        marcar(df[f'flag_{b}'].fillna(0) > 0,                           'flag_galfitm')
        marcar(df[f'e_MAG_{b}'].isna() | (df[f'e_MAG_{b}'] <= 0),       'sin_error')
        marcar((df[f'MAG_{b}'] < MAG_MIN) | (df[f'MAG_{b}'] > MAG_MAX), 'mag_fuera')
        marcar((df[f'RE_{b}'] < RE_MIN) | (df[f'RE_{b}'] > RE_MAX),     're_fuera')
        marcar((df[f'N_{b}'] < N_MIN) | (df[f'N_{b}'] > N_MAX),         'n_limite')
        marcar(df[f'AR_{b}'] < AR_MIN,                                  'ar_bajo')

        col = df[f'MAG_{b}'] - med_gal
        c0  = col[tiene & (malo == '')].median()
        sig = 1.4826 * (col[tiene & (malo == '')] - c0).abs().median()
        lim = max(COLOR_NSIGMA * sig, COLOR_MIN_DESV)
        marcar((col - c0).abs() > lim, 'color_outlier')

        if f'frac_vacia_{b}' in df:
            marcar(df[f'frac_vacia_{b}'] > FRAC_VACIA_MAX, 'imagen_vacia')
            marcar(df[f'fondo_{b}'] > FONDO_MAX,           'fondo_malo')

        df_clean[f'banda_ok_{b}'] = tiene & (malo == '')
        for p in PARAMS + ['SKY']:
            for c in (f'{p}_{b}', f'e_{p}_{b}'):
                if c in df_clean:
                    df_clean.loc[tiene & (malo != ''), c] = np.nan
        for idx in df.index[malo != '']:
            motivos_banda[df.at[idx, 'objectId']].append(f'{b}:{malo[idx]}')

        print(f'  banda {b}: {tiene.sum():5d} con datos, {(malo != "").sum():4d} descartadas '
              f'{ {k: int(v) for k, v in malo[malo != ""].value_counts().items()} }')

    df_clean['bandas_malas'] = df_clean.objectId.map(lambda o: ';'.join(motivos_banda[o]))
    cols_ok = [c for c in df_clean if c.startswith('banda_ok_')]
    df_clean['n_bandas_ok'] = df_clean[cols_ok].sum(axis=1)

    # Una galaxia con bandas malas: las demas bandas tambien pueden estar afectadas,
    # porque GalfitM ajusta todas juntas (polinomio en lambda).
    df_clean['revisar'] = (~df_clean.galaxia_ok) | (df_clean.bandas_malas != '')
    df_clean['RE_arcsec_i'] = df_clean.get('RE_i') * PIXSCALE

    df_clean.to_csv(os.path.join(OUT_DIR, 'galfitm_limpio.csv'), index=False)
    rev = df_clean[df_clean.revisar][['objectId', 'chi2nu', 'bandas',
                                      'motivo_galaxia', 'bandas_malas']]
    rev.to_csv(os.path.join(OUT_DIR, 'galaxias_a_revisar.csv'), index=False)

    buenas = df_clean[df_clean.galaxia_ok]
    print()
    print('=' * 60)
    print(f'  Galaxias OK (pasan cortes de galaxia): {len(buenas)} / {len(df)}')
    print(f'  Descartadas: { {k: int(v) for k, v in df.motivo_galaxia[~df.galaxia_ok].value_counts().items()} }')
    print(f'  OK pero con alguna banda mala        : {(buenas.bandas_malas != "").sum()}')
    print(f'  Totalmente limpias                   : {(buenas.bandas_malas == "").sum()}')
    print('=' * 60)


    # =============================================================================
    # 3. PLOTS
    # =============================================================================
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.3})
    bs = [b for b in BANDAS if f'MAG_{b}' in buenas]
    lam = np.array([LAMBDA_EFF[b] for b in bs])

    def guardar(fig, nombre):
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, nombre), dpi=120)
        plt.close(fig)

    # --- 01 chi2
    fig, ax = plt.subplots(figsize=(7, 4))
    x = df.chi2nu[df.estado == 'ok'].clip(1e-2, 1e2)
    ax.hist(x, bins=np.logspace(-2, 2, 80), color='gray')
    for v in (CHI2_MIN, CHI2_MAX):
        ax.axvline(v, color='r', ls='--')
    ax.set_xscale('log'); ax.set_xlabel(r'$\chi^2_\nu$'); ax.set_ylabel('Number of galaxies')
    ax.set_title(f'chi2nu (cuts at {CHI2_MIN} and {CHI2_MAX})')
    guardar(fig, '01_chi2.png')

    # --- 02 histogramas por banda: antes (gris) vs despues (color)
    specs = [('MAG', np.linspace(12, 26, 60), 'Magnitude', False),
             ('RE',  np.logspace(-0.5, 3, 60), 'R_e [pix]', True),
             ('N',   np.linspace(0, 10, 60), 'Sersic index n', False),
             ('AR',  np.linspace(0, 1, 50), 'b/a', False)]
    fig, axes = plt.subplots(len(specs), len(bs), figsize=(3 * len(bs), 10), sharex='row')
    for r, (p, bins, lab, logx) in enumerate(specs):
        for c, b in enumerate(bs):
            ax = axes[r, c]
            ax.hist(df[f'{p}_{b}'].dropna(), bins=bins, color='lightgray', label='all')
            ax.hist(buenas[f'{p}_{b}'].dropna(), bins=bins, color=COLORES[b],
                    alpha=0.8, label='after cuts')
            if logx: ax.set_xscale('log')
            if r == 0: ax.set_title(b)
            if c == 0: ax.set_ylabel(lab)
    axes[0, 0].legend(fontsize=7)
    guardar(fig, '02_histogramas_por_banda.png')

    # --- 03 dependencia con lambda (lo central de GalfitM)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.3))
    ref = 'i' if 'i' in bs else bs[len(bs) // 2]
    for ax, p, lab in [(axes[0], 'RE', f'R_e / R_e({ref})'),
                       (axes[1], 'N',  f'n / n({ref})'),
                       (axes[2], 'AR', f'(b/a) / (b/a)({ref})')]:
        rel = np.array([buenas[f'{p}_{b}'] / buenas[f'{p}_{ref}'] for b in bs])
        muestra = np.random.default_rng(1).choice(rel.shape[1], min(150, rel.shape[1]), replace=False)
        for j in muestra:
            ax.plot(lam, rel[:, j], color='gray', alpha=0.08, lw=0.8)
        p16, p50, p84 = np.nanpercentile(rel, [16, 50, 84], axis=1)
        ax.fill_between(lam, p16, p84, color='C0', alpha=0.25, label='16-84%')
        ax.plot(lam, p50, 'o-', color='C0', lw=2, label='median')
        ax.axhline(1, color='k', lw=0.8)
        ax.set_ylim(0.5, 1.5); ax.set_xlabel(r'$\lambda$ [$\AA$]'); ax.set_ylabel(lab)
        for b, l in zip(bs, lam):
            ax.text(l, 1.45, b, ha='center')
    axes[0].legend(fontsize=8)
    fig.suptitle('Wavelength dependence of structural parameters')
    guardar(fig, '03_dependencia_lambda.png')

    # --- 04 colores
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    if {'g', 'r', 'i'} <= set(bs):
        gr = buenas.MAG_g - buenas.MAG_r
        ri = buenas.MAG_r - buenas.MAG_i
        axes[0].scatter(buenas.MAG_r, gr, s=3, alpha=0.4, c='k')
        axes[0].set_xlabel('r'); axes[0].set_ylabel('g - r'); axes[0].set_ylim(-0.5, 2)
        axes[0].set_title('Color-magnitude')
        axes[1].scatter(ri, gr, s=3, alpha=0.4, c='k')
        axes[1].set_xlabel('r - i'); axes[1].set_ylabel('g - r')
        axes[1].set_xlim(-0.5, 1.5); axes[1].set_ylim(-0.5, 2); axes[1].set_title('Color-color')
        sc = axes[2].scatter(gr, buenas.N_r, s=3, alpha=0.5, c=np.log10(buenas.RE_r), cmap='viridis')
        axes[2].set_xlabel('g - r'); axes[2].set_ylabel('n (r)'); axes[2].set_yscale('log')
        axes[2].set_xlim(-0.5, 2); axes[2].set_title('Color vs Sersic n')
        plt.colorbar(sc, ax=axes[2], label='log R_e(r) [pix]')
    guardar(fig, '04_colores.png')

    # --- 05 n vs Re y b/a
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(buenas[f'RE_{ref}'] * PIXSCALE, buenas[f'N_{ref}'], s=3, alpha=0.4, c='k')
    axes[0].set_xscale('log'); axes[0].set_yscale('log')
    axes[0].set_xlabel(f'R_e({ref}) [arcsec]'); axes[0].set_ylabel(f'n({ref})')
    axes[1].scatter(buenas[f'N_{ref}'], buenas[f'AR_{ref}'], s=3, alpha=0.4, c='k')
    axes[1].set_xscale('log'); axes[1].set_xlabel(f'n({ref})'); axes[1].set_ylabel(f'b/a({ref})')
    guardar(fig, '05_n_re_ar.png')

    # --- 06 errores vs magnitud
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for b in bs:
        m = buenas[f'MAG_{b}']
        axes[0].scatter(m, buenas[f'e_MAG_{b}'], s=2, alpha=0.3, color=COLORES[b], label=b)
        axes[1].scatter(m, buenas[f'e_RE_{b}'] / buenas[f'RE_{b}'], s=2, alpha=0.3, color=COLORES[b])
        axes[2].scatter(m, buenas[f'e_N_{b}'] / buenas[f'N_{b}'], s=2, alpha=0.3, color=COLORES[b])
    for ax, lab in zip(axes, ['err MAG', 'err R_e / R_e', 'err n / n']):
        ax.set_yscale('log'); ax.set_xlabel('Magnitude'); ax.set_ylabel(lab)
    axes[0].legend(markerscale=5, fontsize=8)
    fig.suptitle('GalfitM formal errors (usually underestimated)')
    guardar(fig, '06_errores.png')

    # --- 07 calidad de imagen (para calibrar los cortes)
    if any(c.startswith('fondo_') for c in df):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        for b in bs:
            if f'fondo_{b}' in df:
                axes[0].hist(df[f'fondo_{b}'].replace(np.inf, 1e3).clip(1e-2, 1e3).dropna(),
                             bins=np.logspace(-2, 3, 60), histtype='step', color=COLORES[b], label=b)
                axes[1].hist(df[f'frac_vacia_{b}'].dropna(), bins=np.linspace(0, 1, 41),
                             histtype='step', color=COLORES[b], label=b)
        axes[0].axvline(FONDO_MAX, color='r', ls='--'); axes[0].set_xscale('log')
        axes[0].set_xlabel('large-scale background variation / noise'); axes[0].set_yscale('log')
        axes[1].axvline(FRAC_VACIA_MAX, color='r', ls='--')
        axes[1].set_xlabel('fraction of empty pixels'); axes[1].set_yscale('log')
        axes[0].legend()
        guardar(fig, '07_calidad_imagenes.png')

    print(f'\nResultados en: {OUT_DIR}/')



# #############################################################################
# PASO 3: CROSSMATCH CON LSST
# #############################################################################
def paso3():

    # =============================================================================
    # CONFIG
    # =============================================================================
    BASE     = os.path.dirname(os.path.abspath(__file__))
    GALFITM  = os.path.join(BASE, 'analisis_galfitm', 'galfitm_limpio.csv')
    LSST     = LSST_PARQUET   # (configurado arriba)
    OUT_DIR  = os.path.join(BASE, 'analisis_galfitm', 'crossmatch')

    SOLO_GALAXIAS_OK = True     # usa solo las que pasaron los cortes de galaxia
    PIXSCALE  = 0.2             # arcsec / pix
    ZP_NJY    = 31.4            # AB para flujos en nJy
    MAG_LIM   = 40.0            # cModelMag ~80 son valores centinela
    REF       = 'i'             # banda de GalfitM para comparar estructura
                                # (MultiProFit da un solo n para todas las bandas)

    BANDAS     = ['u', 'g', 'r', 'i', 'z', 'y']
    COLORES    = {'u': '#7b2cbf', 'g': '#2a9d8f', 'r': '#e76f51',
                  'i': '#d62828', 'z': '#6d3a1f', 'y': '#333333'}

    os.makedirs(OUT_DIR, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.3})


    def guardar(fig, nombre):
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, nombre), dpi=120)
        plt.close(fig)


    def stats_robustas(d):
        """mediana, sigma robusta (MAD) y fraccion de outliers > 0.5 mag"""
        d = d[np.isfinite(d)]
        if len(d) == 0:
            return np.nan, np.nan, np.nan, 0
        med = np.median(d)
        sig = 1.4826 * np.median(np.abs(d - med))
        return med, sig, np.mean(np.abs(d - med) > 0.5), len(d)


    def mediana_en_bins(x, y, bins):
        idx = np.digitize(x, bins)
        xc, med, p16, p84 = [], [], [], []
        for k in range(1, len(bins)):
            s = y[(idx == k) & np.isfinite(y)]
            if len(s) >= 10:
                xc.append(0.5 * (bins[k - 1] + bins[k]))
                med.append(np.median(s))
                p16.append(np.percentile(s, 16))
                p84.append(np.percentile(s, 84))
        return map(np.array, (xc, med, p16, p84))


    # =============================================================================
    # 1. LEER Y CRUZAR
    # =============================================================================
    gf = pd.read_csv(GALFITM, dtype={'objectId': str}, low_memory=False)
    gf['objectId'] = gf['objectId'].astype('int64')
    if SOLO_GALAXIAS_OK and 'galaxia_ok' in gf:
        gf = gf[gf.galaxia_ok.astype(bool)]
    print(f'GalfitM: {len(gf)} galaxias')

    cols_lsst = (['objectId', 'ID', 'z', 'f_z', 'class', 'source', '_dist_arcsec',
                  'coord_ra', 'coord_dec', 'ebv',
                  'deblend_blendNChild', 'detect_fromBlend', 'detect_isIsolated',
                  'sersic_index', 'sersic_indexErr', 'sersic_reff_major',
                  'sersic_reff_minor', 'sersic_theta', 'sersic_chi2_reduced',
                  'sersic_unknown_flag', 'sersic_no_data_flag', 'shape_flag']
                 + [f'{b}_{c}' for b in BANDAS
                    for c in ('cModelMag', 'cModelMagErr', 'sersicFlux', 'sersicFluxErr')])
    lsst = pd.read_parquet(LSST)
    faltan = [c for c in cols_lsst if c not in lsst]
    if faltan:
        print(f'Aviso: columnas que no estan en el parquet: {faltan}')
    lsst = lsst[[c for c in cols_lsst if c in lsst]].copy()
    lsst['objectId'] = lsst['objectId'].astype('int64')
    print(f'LSST:    {len(lsst)} filas, {lsst.objectId.nunique()} objectId distintos')

    # Un mismo objeto LSST puede tener varios espectros: nos quedamos con uno,
    # prefiriendo z confiable (f_z=1) y el cruce mas cercano.
    n_dup = lsst.objectId.duplicated().sum()
    if n_dup:
        orden = [c for c in ('f_z', '_dist_arcsec') if c in lsst]
        asc   = [False, True][:len(orden)]
        lsst  = lsst.sort_values(orden, ascending=asc).drop_duplicates('objectId')
        print(f'  {n_dup} filas repetidas por objectId -> se deja 1 por objeto')

    df = gf.merge(lsst, on='objectId', how='left', indicator=True)
    n_match = (df._merge == 'both').sum()
    print(f'Cruzadas: {n_match} / {len(gf)}')
    if n_match < len(gf):
        sin = df.loc[df._merge != 'both', 'objectId']
        sin.to_csv(os.path.join(OUT_DIR, 'sin_cruce.csv'), index=False)
        print(f'  {len(sin)} sin cruce -> sin_cruce.csv')
    df = df[df._merge == 'both'].drop(columns='_merge').copy()

    # Magnitudes LSST limpias y magnitud del Sersic de MultiProFit
    for b in BANDAS:
        if f'{b}_cModelMag' in df:
            m = df[f'{b}_cModelMag'].astype(float)
            df[f'{b}_cModelMag'] = m.where((m > 0) & (m < MAG_LIM))
        if f'{b}_sersicFlux' in df:
            f = df[f'{b}_sersicFlux'].astype(float)
            df[f'{b}_sersicMag'] = np.where(f > 0, -2.5 * np.log10(f.where(f > 0)) + ZP_NJY, np.nan)
        if f'MAG_{b}' in df:
            df[f'dmag_cmodel_{b}'] = df[f'MAG_{b}'] - df.get(f'{b}_cModelMag')
            df[f'dmag_sersic_{b}'] = df[f'MAG_{b}'] - df.get(f'{b}_sersicMag')

    # Estructura LSST en las mismas unidades que GalfitM
    if 'sersic_reff_major' in df:
        df['lsst_RE_pix'] = df.sersic_reff_major / PIXSCALE
        df['lsst_AR']     = df.sersic_reff_minor / df.sersic_reff_major

    # Ajuste de LSST confiable (para los plots de estructura)
    lsst_ok = pd.Series(True, index=df.index)
    for c in ('sersic_unknown_flag', 'sersic_no_data_flag', 'shape_flag'):
        if c in df:
            lsst_ok &= ~df[c].fillna(True).astype(bool)
    df['lsst_sersic_ok'] = lsst_ok

    df['objectId'] = df['objectId'].astype(str)    # evita redondeo al abrir en Excel
    df.to_csv(os.path.join(OUT_DIR, 'galfitm_x_lsst.csv'), index=False)


    # =============================================================================
    # 2. RESUMEN NUMERICO DE MAGNITUDES
    # =============================================================================
    filas = []
    for b in BANDAS:
        for tipo, lab in (('cmodel', 'cModel'), ('sersic', 'Sersic LSST')):
            c = f'dmag_{tipo}_{b}'
            if c in df:
                med, sig, fout, n = stats_robustas(df[c].values.astype(float))
                filas.append({'banda': b, 'comparacion': f'GalfitM - {lab}', 'N': n,
                              'offset_mediano': med, 'sigma_robusta': sig,
                              'frac_|d|>0.5': fout})
    res = pd.DataFrame(filas)
    res.to_csv(os.path.join(OUT_DIR, 'resumen_magnitudes.csv'), index=False)
    print('\n' + res.round(3).to_string(index=False))


    # =============================================================================
    # 3. PLOTS
    # =============================================================================
    bs = [b for b in BANDAS if f'dmag_cmodel_{b}' in df and df[f'dmag_cmodel_{b}'].notna().sum() > 10]
    bins_mag = np.arange(12, 24, 0.5)

    def panel_mag(tipo, col_lsst, titulo, nombre):
        fig, axes = plt.subplots(2, len(bs), figsize=(3.3 * len(bs), 6.5),
                                 squeeze=False, sharey='row')
        for c, b in enumerate(bs):
            x = df[col_lsst.format(b=b)].astype(float).values
            y = df[f'MAG_{b}'].astype(float).values
            d = y - x
            ok = np.isfinite(x) & np.isfinite(y)
            med, sig, fout, n = stats_robustas(d[ok])

            ax = axes[0, c]
            ax.scatter(x[ok], y[ok], s=2, alpha=0.3, color=COLORES[b])
            lim = [11, 24]
            ax.plot(lim, lim, 'k-', lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim)
            ax.set_title(f'{b}   (N={n})')
            ax.set_xlabel(f'{b} LSST'); 
            if c == 0: ax.set_ylabel('GalfitM')

            ax = axes[1, c]
            ax.scatter(x[ok], d[ok], s=2, alpha=0.3, color=COLORES[b])
            xc, m50, p16, p84 = mediana_en_bins(x[ok], d[ok], bins_mag)
            ax.plot(xc, m50, 'k-', lw=1.5)
            ax.fill_between(xc, p16, p84, color='k', alpha=0.15)
            ax.axhline(0, color='r', lw=0.8)
            ax.set_ylim(-1.5, 1.5); ax.set_xlim(lim)
            ax.text(0.03, 0.05, f'med={med:+.3f}\nσ={sig:.3f}', transform=ax.transAxes,
                    fontsize=8, va='bottom')
            ax.set_xlabel(f'{b} LSST')
            if c == 0: ax.set_ylabel('GalfitM − LSST [mag]')
        fig.suptitle(titulo)
        guardar(fig, nombre)

    panel_mag('cmodel', '{b}_cModelMag', 'GalfitM vs LSST cModelMag', '08_mag_vs_cmodel.png')
    if any(f'{b}_sersicMag' in df for b in bs):
        panel_mag('sersic', '{b}_sersicMag', 'GalfitM vs LSST Sersic (MultiProFit)',
                  '09_mag_vs_sersic_lsst.png')

    # --- Estructura: R_e, n, b/a
    s = df[df.lsst_sersic_ok]
    if 'lsst_RE_pix' in df and f'RE_{REF}' in df:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
        pares = [(f'RE_{REF}', 'lsst_RE_pix', f'R_e GalfitM ({REF}) [arcsec]', 'R_e (major axis) LSST [arcsec]', True, PIXSCALE),
                 (f'N_{REF}',  'sersic_index', f'n GalfitM ({REF})',           'n LSST',                   True, 1.0),
                 (f'AR_{REF}', 'lsst_AR',      f'b/a GalfitM ({REF})',          'b/a LSST',                 False, 1.0)]
        for ax, (cg, cl, lg, ll, logs, esc) in zip(axes, pares):
            xg = s[cg].astype(float) * esc
            xl = s[cl].astype(float) * esc
            ok = np.isfinite(xg) & np.isfinite(xl) & (xg > 0) & (xl > 0)
            ax.scatter(xl[ok], xg[ok], s=3, alpha=0.3, c='k')
            lo, hi = np.nanpercentile(np.r_[xl[ok], xg[ok]], [0.5, 99.5])
            ax.plot([lo, hi], [lo, hi], 'r-', lw=1)
            if logs:
                ax.set_xscale('log'); ax.set_yscale('log')
                r = np.log10(xg[ok] / xl[ok])
                txt = f'median log(GF/LSST)={np.median(r):+.3f}\nσ={1.4826*np.median(np.abs(r-np.median(r))):.3f} dex'
            else:
                r = xg[ok] - xl[ok]
                txt = f'median Δ={np.median(r):+.3f}\nσ={1.4826*np.median(np.abs(r-np.median(r))):.3f}'
            ax.text(0.03, 0.97, f'N={ok.sum()}\n' + txt, transform=ax.transAxes, va='top', fontsize=8)
            ax.set_xlabel(ll); ax.set_ylabel(lg)
        # limites del n de MultiProFit
        axes[1].axvline(0.5, color='gray', ls=':'); axes[1].axvline(6, color='gray', ls=':')
        fig.suptitle(f'Structure: GalfitM ({REF} band) vs MultiProFit (LSST, unflagged)')
        guardar(fig, '10_estructura_vs_lsst.png')

    # --- Diferencia de magnitud vs n: las Sersic de n alto tienen alas extensas,
    #     asi que la magnitud total suele ser mas brillante que cModel
    if f'dmag_cmodel_{REF}' in df:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        for ax, (tipo, lab) in zip(axes, (('cmodel', 'cModel'), ('sersic', 'Sersic LSST'))):
            c = f'dmag_{tipo}_{REF}'
            if c not in df:
                continue
            x = df[f'N_{REF}'].astype(float).values
            y = df[c].astype(float).values
            ok = np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[ok], y[ok], s=2, alpha=0.3, c='k')
            xc, m50, p16, p84 = mediana_en_bins(x[ok], y[ok], np.arange(0.2, 8.5, 0.5))
            ax.plot(xc, m50, 'r-', lw=2, label='median')
            ax.fill_between(xc, p16, p84, color='r', alpha=0.2)
            ax.axhline(0, color='k', lw=0.8)
            ax.set_ylim(-1.5, 1.0)
            ax.set_xlabel(f'n GalfitM ({REF})'); ax.set_ylabel(f'GalfitM − {lab} ({REF}) [mag]')
            ax.legend(fontsize=8)
        fig.suptitle('Magnitude difference vs Sersic index')
        guardar(fig, '11_dmag_vs_n.png')

    print(f'\nResultados en: {OUT_DIR}/')



# #############################################################################
# PASO 4: MUESTRA CONFIABLE (A/B/C/D)
# #############################################################################
def paso4():

    # =============================================================================
    # CONFIG
    # =============================================================================
    BASE        = os.path.dirname(os.path.abspath(__file__))
    CSV_IN      = os.path.join(BASE, 'analisis_galfitm', 'crossmatch', 'galfitm_x_lsst.csv')
    GALFITM_DIR = os.path.join(BASE, 'galfitm_inputs')
    OUT_DIR     = os.path.join(BASE, 'analisis_galfitm', 'muestra')
    METRICAS    = os.path.join(OUT_DIR, 'metricas_imagen.csv')
    RECALCULAR  = RECALCULAR_METRICAS   # (configurado arriba)
    PIXSCALE    = 0.2
    ZP          = 31.4           # zeropoint de los feedme (nJy -> AB)

    BANDAS_REF    = ['i', 'r', 'z', 'g', 'y']
    MIN_BANDAS_OK = 3
    BANDAS_OBLIGATORIAS = ['g', 'r', 'i']   # bandas que TODAS las galaxias deben tener buenas
                                       # g, r: clasificacion (g-r, n_r);  i: masa de Taylor (g-i, M_i)

    # Espectroscopia
    F_Z_OK    = [1]
    CLASES_OK = ['GALAXY']
    DIST_MAX  = 1.0

    # Acuerdo GalfitM-LSST (respecto al offset mediano, medido de los datos)
    # Cortes SUAVES: LSST no se toma como verdad; solo se marca "desacuerdo" cuando
    # la diferencia es enorme. Tolerancia = max(minimo absoluto, N_SIGMA_LSST * sigma
    # robusta de las diferencias), siempre medida respecto al offset mediano.
    # Aun asi, una galaxia en desacuerdo NO se elimina si GalfitM reproduce la
    # imagen (categoria B); solo va a C si ademas falla el arbitro de la imagen.
    N_SIGMA_LSST = 5.0
    TOL_MAG = 1.0     # mag   (un factor 2.5 en flujo)
    TOL_RE  = 0.30    # dex   (un factor 2 en R_e)
    TOL_N   = 0.30    # dex   (un factor 2 en n)
    TOL_AR  = 0.25
    N_LIM_LSST = (0.5, 6.0)

    # Arbitro: la imagen
    R_APER    = 2.0   # apertura eliptica de radio R_APER * R_e (forma de GalfitM)
    RFF_MAX   = 0.10  # residual flux fraction maxima
    FLUJO_TOL = 0.10  # |flujo modelo / flujo dato - 1| maximo en la apertura
    R_TOTAL   = 5.0   # curva de crecimiento: flujo "total" dentro de R_TOTAL * R_e(GalfitM)

    # Problemas internos de LSST (solo informativos)
    NCHILD_FRAGMENTADA = 10
    CHI2_LSST_MAX      = 5.0

    os.makedirs(OUT_DIR, exist_ok=True)
    plt.rcParams.update({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.3})


    def mad(x):
        x = np.asarray(x, float); x = x[np.isfinite(x)]
        return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


    # =============================================================================
    # 1. LEER Y ELEGIR BANDA DE REFERENCIA
    # =============================================================================
    df = pd.read_csv(CSV_IN, dtype={'objectId': str}, low_memory=False).copy()
    print(f'Galaxias en la tabla cruzada: {len(df)}')
    # columnas opcionales de LSST (solo informativas): si faltan, valores neutros
    for c, v in (('shape_flag', False), ('deblend_blendNChild', np.nan),
                 ('sersic_chi2_reduced', np.nan), ('z', np.nan)):
        if c not in df:
            df[c] = v

    ref = pd.Series(None, index=df.index, dtype=object)
    for b in reversed(BANDAS_REF):            # la primera de la lista gana
        if f'MAG_{b}' in df:
            ref[df[f'MAG_{b}'].notna()] = b
    df['banda_ref'] = ref

    def col_ref(fmt):
        out = np.full(len(df), np.nan)
        for b in BANDAS_REF:
            c = fmt.format(b=b)
            if c in df:
                m = (df.banda_ref == b).values
                out[m] = pd.to_numeric(df.loc[m, c], errors='coerce').values
        return out

    for p in ('MAG', 'RE', 'N', 'AR', 'PA', 'XC', 'YC'):
        df[f'{p}_ref'] = col_ref(p + '_{b}')
    df['lsst_mag_ref'] = col_ref('{b}_sersicMag')


    # =============================================================================
    # 2. ARBITRO: METRICAS MEDIDAS DIRECTAMENTE EN LA IMAGEN
    # =============================================================================
    def radio_eliptico(shape, xc, yc, re, q, pa):
        """Radio eliptico en unidades de R_e, convencion de GALFIT:
           PA en grados desde +y hacia -x; xc, yc en pixeles base 1."""
        yy, xx = np.indices(shape)
        dx, dy = xx + 1 - xc, yy + 1 - yc
        t = np.radians(pa)
        a = -dx * np.sin(t) + dy * np.cos(t)      # eje mayor
        b =  dx * np.cos(t) + dy * np.sin(t)      # eje menor
        return np.sqrt(a**2 + (b / max(q, 0.05))**2) / max(re, 0.5)


    def medir(fila):
        from astropy.io import fits
        oid, b = fila.objectId, fila.banda_ref
        carpeta = dir_galaxia(oid, getattr(fila, 'ruta_carpeta', None))
        out = {'objectId': oid}
        with fits.open(os.path.join(carpeta, f'{oid}_galfitm_out.fits')) as h:
            nom = {x.name.upper(): x for x in h}
            dato = nom[f'INPUT_{b.upper()}'].data.astype(float)
            mod  = nom[f'MODEL_{b.upper()}'].data.astype(float)
            res  = nom[f'RESIDUAL_{b.upper()}'].data.astype(float)

        malo = ~np.isfinite(dato) | (dato == 0)
        pm = os.path.join(carpeta, f'{oid}_mask.fits')
        if os.path.exists(pm):
            m = fits.getdata(pm)
            if m.shape == dato.shape:
                malo |= (m != 0)

        sig = None
        ps = os.path.join(carpeta, f'{oid}_sig_{b}.fits')
        if os.path.exists(ps):
            s = fits.getdata(ps).astype(float)
            if s.shape == dato.shape:
                sig = s

        r = radio_eliptico(dato.shape, fila.XC_ref, fila.YC_ref,
                           fila.RE_ref, fila.AR_ref, fila.PA_ref)
        ap = (r <= R_APER) & ~malo
        out['frac_apertura_enmascarada'] = 1 - ap.sum() / max((r <= R_APER).sum(), 1)

        # Cielo: mediana en un anillo lejano. Se resta igual al dato y al modelo
        anillo = (r > 4 * R_APER) & ~malo
        cielo = np.median(dato[anillo]) if anillo.sum() > 100 else 0.0

        f_dato = np.sum(dato[ap] - cielo)
        f_mod  = np.sum(mod[ap] - cielo)
        out['ratio_flujo_ap'] = f_mod / f_dato if f_dato > 0 else np.nan
        out['mag_aper_dato']  = -2.5 * np.log10(f_dato) + ZP if f_dato > 0 else np.nan

        # Residual Flux Fraction (Hoyos et al. 2011): residuo en exceso del ruido
        ruido = 0.8 * np.sum(sig[ap]) if sig is not None else 0.0
        out['rff'] = (np.sum(np.abs(res[ap])) - ruido) / f_dato if f_dato > 0 else np.nan

        # --- Curva de crecimiento en la IMAGEN (sin modelo) para juzgar R_e ---
        # Pixeles enmascarados (vecinos) se rellenan con el modelo de GalfitM.
        # Misma elipse (q, PA de GalfitM); solo cambia el radio que se prueba.
        lleno = np.where(malo, mod, dato) - cielo
        dentro = r <= R_TOTAL
        rr, ff = r[dentro], lleno[dentro]
        orden = np.argsort(rr)
        rr, cum = rr[orden], np.cumsum(ff[orden])
        f_tot = cum[-1] if len(cum) else np.nan
        if f_tot > 0:
            frac = lambda R: np.interp(R, rr, cum) / f_tot
            out['f_dentro_re_galfitm'] = frac(1.0)
            re_lsst = getattr(fila, 'sersic_reff_major', np.nan)
            if pd.notna(re_lsst) and re_lsst > 0:
                out['f_dentro_re_lsst'] = frac(re_lsst / PIXSCALE / fila.RE_ref)
            k = np.searchsorted(cum, 0.5 * f_tot)
            out['r50_imagen_pix'] = rr[min(k, len(rr) - 1)] * fila.RE_ref
        # fraccion esperada para un perfil de Sersic con el n de GalfitM
        from scipy.special import gammainc, gammaincinv
        n_ = float(fila.N_ref)
        if np.isfinite(n_) and n_ > 0:
            b_ = gammaincinv(2 * n_, 0.5)
            out['f_esperada_sersic'] = 0.5 / gammainc(2 * n_, b_ * R_TOTAL**(1 / n_))
        return out


    candidatas = df[df.galaxia_ok.astype(bool) & df.banda_ref.notna()]
    # Cache incremental: solo se miden las galaxias que faltan
    met = pd.DataFrame({'objectId': pd.Series([], dtype=str)})
    if os.path.exists(METRICAS) and not RECALCULAR:
        met = pd.read_csv(METRICAS, dtype={'objectId': str})
        if 'f_dentro_re_galfitm' not in met:
            print('Las metricas guardadas son de una version anterior: se recalculan todas')
            met = pd.DataFrame({'objectId': pd.Series([], dtype=str)})
    faltan = candidatas[~candidatas.objectId.isin(met.objectId)]
    if len(faltan) and hay_imagenes(df):
        print(f'Metricas de imagen: {len(met)} guardadas, midiendo {len(faltan)} nuevas...')
        filas = []
        for k, (_, f) in enumerate(faltan.iterrows(), 1):
            try:
                filas.append(medir(f))
            except Exception as e:
                filas.append({'objectId': f.objectId, 'error_imagen': str(e)[:80]})
            if k % 500 == 0:
                print(f'  {k}/{len(faltan)}')
        met = pd.concat([met, pd.DataFrame(filas)], ignore_index=True)
        met.to_csv(METRICAS, index=False)
    elif len(faltan):
        print('AVISO: no encuentro las carpetas galfitm_inputs; sin arbitro de imagen '
              'esas galaxias quedan en C si no coinciden con LSST')
    else:
        print(f'Metricas de imagen leidas de {METRICAS} (todas ya medidas)')
    df = df.merge(met, on='objectId', how='left')
    for c in ('rff', 'ratio_flujo_ap', 'mag_aper_dato', 'f_dentro_re_galfitm',
              'f_dentro_re_lsst', 'r50_imagen_pix', 'f_esperada_sersic'):
        if c not in df:
            df[c] = np.nan


    # =============================================================================
    # 3. COMPARACION GALFITM - LSST (dos medidas, ninguna es la verdad)
    # =============================================================================
    num = lambda c: pd.to_numeric(df[c], errors='coerce')
    df['d_mag']   = df.MAG_ref - df.lsst_mag_ref
    df['d_logre'] = np.log10(df.RE_ref * PIXSCALE / num('sersic_reff_major'))
    df['d_logn']  = np.log10(df.N_ref / num('sersic_index'))
    df['d_ar']    = df.AR_ref - num('lsst_AR')

    # Tercera opinion sobre el flujo, sin modelo: cada uno contra la apertura
    df['d_mag_gf_vs_aper']   = df.MAG_ref - df.mag_aper_dato
    df['d_mag_lsst_vs_aper'] = df.lsst_mag_ref - df.mag_aper_dato

    comparable = (df.lsst_sersic_ok.astype(bool) & df.banda_ref.notna()
                  & df.lsst_mag_ref.notna())

    offs, sigs = {}, {}
    for c in ('d_mag', 'd_logre', 'd_logn', 'd_ar'):
        offs[c] = np.nanmedian(df.loc[comparable, c]); sigs[c] = mad(df.loc[comparable, c])
    print('\nDiferencias sistematicas GalfitM - LSST (mediana, sigma robusta):')
    for c in offs:
        print(f'  {c:8s} {offs[c]:+.3f}   sigma={sigs[c]:.3f}')

    n_lim = ((num('sersic_index') <= N_LIM_LSST[0] + 0.01) |
             (num('sersic_index') >= N_LIM_LSST[1] - 0.01))
    tol = {c: max(t, N_SIGMA_LSST * sigs[c]) for c, t in
           (('d_mag', TOL_MAG), ('d_logre', TOL_RE), ('d_logn', TOL_N), ('d_ar', TOL_AR))}
    print('Tolerancias usadas (desacuerdo solo si se supera): ' +
          ', '.join(f'{c} {v:.2f}' for c, v in tol.items()))
    ok_mag = (df.d_mag   - offs['d_mag']).abs()   < tol['d_mag']
    ok_re  = (df.d_logre - offs['d_logre']).abs() < tol['d_logre']
    ok_n   = n_lim | ((df.d_logn - offs['d_logn']).abs() < tol['d_logn'])
    ok_ar  = (df.d_ar    - offs['d_ar']).abs()    < tol['d_ar']
    acuerdo = comparable & ok_mag & ok_re & ok_n & ok_ar

    df['comparacion_lsst'] = np.where(~comparable, 'no_comparable',
                             np.where(acuerdo, 'acuerdo', 'desacuerdo'))

    # Problemas internos de LSST (informativo)
    shape  = df.shape_flag.fillna(False).astype(bool)
    frag   = num('deblend_blendNChild').fillna(0) >= NCHILD_FRAGMENTADA
    chi2l  = num('sersic_chi2_reduced').fillna(0) > CHI2_LSST_MAX
    df['problemas_lsst'] = [';'.join(t for t, v in (('shape_flag', a), ('fragmentada', b),
                                                    ('n_en_limite', c), ('chi2_alto', d)) if v)
                            for a, b, c, d in zip(shape, frag, n_lim, chi2l)]

    imagen_ok = (df.rff < RFF_MAX) & ((df.ratio_flujo_ap - 1).abs() < FLUJO_TOL)
    df['galfitm_reproduce_imagen'] = np.where(df.rff.isna(), 'sin_medir',
                                     np.where(imagen_ok, 'si', 'no'))


    # =============================================================================
    # 4. CATEGORIAS
    # =============================================================================
    # Banda "buena" = tiene MAG y n (los cortes por banda de analizar_galfitm.py ya
    # pusieron NaN a las bandas malas, asi que esto exige que pasen esos cortes)
    tiene_obligatorias = pd.Series(True, index=df.index)
    for b in BANDAS_OBLIGATORIAS:
        tiene_obligatorias &= df[f'MAG_{b}'].notna() & df[f'N_{b}'].notna()

    c_galfitm = (df.galaxia_ok.astype(bool) & df.banda_ref.notna()
                 & (df.n_bandas_ok >= MIN_BANDAS_OK) & tiene_obligatorias)
    c_spec = df.f_z.isin(F_Z_OK) & df['class'].isin(CLASES_OK) & (df._dist_arcsec <= DIST_MAX)

    cat = np.where(~(c_galfitm & c_spec), 'D',
          np.where(acuerdo, 'A',
          np.where(imagen_ok, 'B', 'C')))
    df['categoria'] = cat

    def motivo(i):
        m = []
        if not c_galfitm.iloc[i]:
            if not bool(df.galaxia_ok.iloc[i]):      m.append(f'galfitm:{df.motivo_galaxia.iloc[i]}')
            else:
                if pd.isna(df.banda_ref.iloc[i]):    m.append('sin_banda_ref')
                if df.n_bandas_ok.iloc[i] < MIN_BANDAS_OK: m.append('pocas_bandas_ok')
                for b in BANDAS_OBLIGATORIAS:
                    if pd.isna(df[f'MAG_{b}'].iloc[i]) or pd.isna(df[f'N_{b}'].iloc[i]):
                        m.append(f'falta_banda_{b}')
        if not c_spec.iloc[i]:
            if df.f_z.iloc[i] not in F_Z_OK:         m.append(f'f_z={df.f_z.iloc[i]}')
            if df['class'].iloc[i] not in CLASES_OK: m.append(f'clase={df["class"].iloc[i]}')
            if df._dist_arcsec.iloc[i] > DIST_MAX:   m.append('dist_grande')
        if cat[i] in ('B', 'C'):
            if comparable.iloc[i]:
                for ok, t in ((ok_mag, 'dif_mag'), (ok_re, 'dif_re'), (ok_n, 'dif_n'), (ok_ar, 'dif_ar')):
                    if not ok.iloc[i]:
                        m.append(t)
            else:
                m.append('lsst_no_comparable')
            if cat[i] == 'C':
                m.append(f'reproduce_imagen={df.galfitm_reproduce_imagen.iloc[i]}')
        return ';'.join(m)
    df['motivo_categoria'] = [motivo(i) for i in range(len(df))]

    df.to_csv(os.path.join(OUT_DIR, 'muestra_clasificada.csv'), index=False)
    df[df.categoria == 'A'].to_csv(os.path.join(OUT_DIR, 'muestra_A.csv'), index=False)
    df[df.categoria.isin(['A', 'B'])].to_csv(os.path.join(OUT_DIR, 'muestra_AB.csv'), index=False)

    print('\n' + '=' * 64)
    desc = {'A': 'coincide con LSST', 'B': 'no coincide, GalfitM reproduce la imagen',
            'C': 'no coincide, GalfitM no reproduce la imagen', 'D': 'descartada'}
    for k in 'ABCD':
        print(f'  {k}  {(df.categoria == k).sum():5d}   {desc[k]}')
    print('=' * 64)
    print(f'  Sin alguna banda obligatoria {BANDAS_OBLIGATORIAS}: '
          f'{(df.galaxia_ok.astype(bool) & ~tiene_obligatorias).sum()} (van a D)')

    # Quien se equivoca en los desacuerdos? Comparacion con la apertura (sin modelo).
    # Un modelo de flujo TOTAL no deberia ser mas debil que el flujo dentro de 2 R_e.
    des = (df.comparacion_lsst == 'desacuerdo') & c_galfitm & c_spec
    s = df[des & df.mag_aper_dato.notna()]
    if len(s) > 10:
        print(f'\nEn {len(s)} desacuerdos con apertura medida:')
        print(f'  LSST mas debil que la apertura (> 0.1 mag)    : '
              f'{(s.d_mag_lsst_vs_aper > 0.1).mean()*100:5.1f}%  -> LSST pierde flujo')
        print(f'  GalfitM mas debil que la apertura (> 0.1 mag) : '
              f'{(s.d_mag_gf_vs_aper > 0.1).mean()*100:5.1f}%  -> GalfitM pierde flujo')

    for k in 'BC':
        cnt = Counter(t for x in df.problemas_lsst[df.categoria == k] for t in x.split(';') if t)
        print(f'Problemas internos de LSST en {k}: {dict(cnt)}')


    # =============================================================================
    # 5. PLOTS
    # =============================================================================
    # ¿Quien acierta el radio efectivo? Fraccion de luz de la IMAGEN dentro de cada R_e
    s = df[(c_galfitm & c_spec) & df.f_dentro_re_galfitm.notna()]
    if len(s) > 10:
        disc = s[s.N_ref < 2]          # perfiles casi exponenciales: el flujo total converge
        print('\nRadio efectivo vs la imagen (galaxias con n < 2, N = %d):' % len(disc))
        print(f'  fraccion esperada dentro de R_e          : {disc.f_esperada_sersic.median():.3f}')
        print(f'  fraccion medida dentro de R_e de GalfitM : {disc.f_dentro_re_galfitm.median():.3f}')
        if disc.f_dentro_re_lsst.notna().any():
            print(f'  fraccion medida dentro de R_e de LSST    : {disc.f_dentro_re_lsst.median():.3f}')
        print(f'  log(r50 imagen / R_e GalfitM): {np.nanmedian(np.log10(disc.r50_imagen_pix / disc.RE_ref)):+.3f} dex')
        re_l = pd.to_numeric(disc.sersic_reff_major, errors='coerce') / PIXSCALE
        print(f'  log(r50 imagen / R_e LSST)   : {np.nanmedian(np.log10(disc.r50_imagen_pix / re_l)):+.3f} dex')
        print('  (r50 incluye el ensanchamiento de la PSF, asi que deberia ser algo mayor que R_e)')

        fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
        bins = np.linspace(0, 1, 51)
        axes[0].hist(disc.f_dentro_re_galfitm.clip(0, 1), bins=bins, histtype='step', lw=2,
                     color='#2e86c1', label='within $R_e$ (GalfitM)')
        axes[0].hist(disc.f_dentro_re_lsst.clip(0, 1).dropna(), bins=bins, histtype='step', lw=2,
                     color='#c0392b', label='within $R_e$ (LSST)')
        axes[0].axvline(disc.f_esperada_sersic.median(), color='k', ls='--',
                        label='expected (Sérsic)')
        axes[0].set_xlabel(f'fraction of image flux (total = {R_TOTAL:g} $R_e$)')
        axes[0].set_title('Galaxies with n < 2', fontsize=10); axes[0].legend(fontsize=8)
        x = s.r50_imagen_pix * PIXSCALE
        axes[1].scatter(x, s.RE_ref * PIXSCALE, s=3, alpha=0.4, color='#2e86c1', label='GalfitM')
        axes[1].scatter(x, pd.to_numeric(s.sersic_reff_major, errors='coerce'), s=3, alpha=0.4,
                        color='#c0392b', label='LSST')
        lim = [np.nanpercentile(x, 1) * 0.7, np.nanpercentile(x, 99) * 1.5]
        axes[1].plot(lim, lim, 'k-', lw=0.8)
        axes[1].set_xscale('log'); axes[1].set_yscale('log')
        axes[1].set_xlabel('half-light radius measured on the image [arcsec]')
        axes[1].set_ylabel('$R_e$ [arcsec]'); axes[1].legend(fontsize=8, markerscale=4)
        axes[2].scatter(s.N_ref, s.f_dentro_re_galfitm, s=3, alpha=0.4, color='#2e86c1', label='GalfitM')
        axes[2].scatter(s.N_ref, s.f_dentro_re_lsst, s=3, alpha=0.4, color='#c0392b', label='LSST')
        nn = np.linspace(0.3, 8, 100)
        from scipy.special import gammainc, gammaincinv
        bb = gammaincinv(2 * nn, 0.5)
        axes[2].plot(nn, 0.5 / gammainc(2 * nn, bb * R_TOTAL**(1 / nn)), 'k--', label='expected')
        axes[2].set_xscale('log'); axes[2].set_ylim(0, 1)
        axes[2].set_xlabel('n (GalfitM)'); axes[2].set_ylabel('flux fraction within $R_e$')
        axes[2].legend(fontsize=8, markerscale=4)
        fig.suptitle('Which effective radius matches the image? (model-independent test)')
        fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, '21_radio_efectivo.png'), dpi=120)
        plt.close(fig)

    COL = {'A': '#d4a017', 'B': '#2e86c1', 'C': '#c0392b', 'D': '#555555'}
    base = c_galfitm & c_spec

    # --- 12 consistencia
    fig, axes = plt.subplots(1, 4, figsize=(17, 4))
    specs = [('d_mag', tol['d_mag'], np.linspace(-2, 1, 90), 'GalfitM − LSST Sérsic [mag]'),
             ('d_logre', tol['d_logre'], np.linspace(-0.4, 1.0, 90), 'log(R_e GF / R_e LSST)'),
             ('d_logn', tol['d_logn'], np.linspace(-0.8, 0.8, 90), 'log(n GF / n LSST)'),
             ('d_ar', tol['d_ar'], np.linspace(-0.5, 0.5, 90), 'b/a GF − b/a LSST')]
    for ax, (c, tol, bins, lab) in zip(axes, specs):
        x = df.loc[base & comparable, c].dropna()
        ax.hist(x, bins=bins, color='gray')
        for v in (offs[c] - tol, offs[c] + tol):
            ax.axvline(v, color='r', ls='--')
        ax.set_title(f'offset={offs[c]:+.3f}  σ={sigs[c]:.3f}', fontsize=9)
        ax.set_xlabel(lab); ax.set_yscale('log')
    fig.suptitle('GalfitM − LSST differences (two independent measurements, neither is ground truth)')
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, '12_consistencia_lsst.png'), dpi=120)
    plt.close(fig)

    # --- 13 arbitro: la imagen decide
    if df.rff.notna().sum() > 10:
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
        for k, col, lab in (('acuerdo', '#d4a017', 'agree'),
                            ('desacuerdo', '#c0392b', 'disagree'),
                            ('no_comparable', '#7f8c8d', 'LSST not comparable')):
            s = df[base & (df.comparacion_lsst == k)]
            axes[0].hist(s.rff.clip(-0.1, 0.5).dropna(), bins=np.linspace(-0.1, 0.5, 50),
                         histtype='step', lw=2, color=col, density=True, label=f'{lab} ({len(s)})')
            axes[1].hist(s.ratio_flujo_ap.clip(0.5, 1.5).dropna(), bins=np.linspace(0.5, 1.5, 50),
                         histtype='step', lw=2, color=col, density=True)
        axes[0].axvline(RFF_MAX, color='k', ls='--')
        axes[0].set_xlabel(f'GalfitM RFF (aperture {R_APER:g} R_e)'); axes[0].legend(fontsize=8)
        for v in (1 - FLUJO_TOL, 1 + FLUJO_TOL):
            axes[1].axvline(v, color='k', ls='--')
        axes[1].set_xlabel(f'model flux / data flux (aperture {R_APER:g} R_e)')
        s = df[base & comparable]
        axes[2].scatter(s.d_mag_lsst_vs_aper, s.d_mag_gf_vs_aper, s=3, alpha=0.3,
                        c=np.where(s.comparacion_lsst == 'acuerdo', '#d4a017', '#c0392b'))
        axes[2].axhline(0, color='k', lw=0.8); axes[2].axvline(0, color='k', lw=0.8)
        axes[2].set_xlim(-2, 1.5); axes[2].set_ylim(-2, 1.5)
        axes[2].set_xlabel('LSST Sérsic − aperture [mag]')
        axes[2].set_ylabel('GalfitM − aperture [mag]')
        axes[2].set_title('A total flux must be ≤ 0 (brighter than the aperture)', fontsize=9)
        fig.suptitle('Model-independent arbiter: the image itself')
        fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, '13_arbitro_imagen.png'), dpi=120)
        plt.close(fig)

    # --- 14 propiedades por categoria
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    props = [('MAG_ref', np.linspace(12, 22, 40), 'GalfitM magnitude', False),
             ('z', np.linspace(0, 0.1, 40), 'z', False),
             ('RE_ref', np.logspace(0.5, 2.5, 40), 'R_e GalfitM [pix]', True),
             ('N_ref', np.linspace(0, 8, 40), 'n GalfitM', False),
             ('AR_ref', np.linspace(0, 1, 30), 'b/a GalfitM', False),
             ('deblend_blendNChild', np.arange(0, 41), 'Number of deblended children (LSST)', False)]
    for ax, (c, bins, lab, logx) in zip(axes.ravel(), props):
        for k in 'ABC':
            x = pd.to_numeric(df.loc[df.categoria == k, c], errors='coerce').dropna()
            if len(x):
                ax.hist(x, bins=bins, histtype='step', lw=2, color=COL[k], density=True,
                        label=f'{k} ({len(x)})')
        if logx:
            ax.set_xscale('log')
        ax.set_xlabel(lab)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Properties by category (normalized)')
    fig.tight_layout(); fig.savefig(os.path.join(OUT_DIR, '14_categorias.png'), dpi=120)
    plt.close(fig)

    print(f'\nResultados en: {OUT_DIR}/')



# =============================================================================
# EJECUCION
# =============================================================================
PASOS = {'2': ('Cortes de calidad de GalfitM', paso2),
         '3': ('Crossmatch con LSST', paso3),
         '4': ('Muestra confiable', paso4)}

REQUISITOS = {'2': os.path.join(DIR_MUESTRA, 'galfitm_resultados.csv'),
              '3': os.path.join(DIR_MUESTRA, 'analisis_galfitm', 'galfitm_limpio.csv'),
              '4': os.path.join(DIR_MUESTRA, 'analisis_galfitm', 'crossmatch', 'galfitm_x_lsst.csv')}

if __name__ == '__main__':
    pedidos = sys.argv[1:] or ['2', '3', '4']
    for p in pedidos:
        if p not in PASOS:
            sys.exit(f'Paso desconocido: {p} (usa 2, 3 o 4)')
    if '3' in pedidos and not os.path.exists(LSST_PARQUET):
        sys.exit(f'No encuentro el parquet de LSST: {LSST_PARQUET}\n'
                 f'Ajusta LSST_PARQUET al inicio del script.')

    for p in pedidos:
        titulo, funcion = PASOS[p]
        if not os.path.exists(REQUISITOS[p]):
            sys.exit(f'\nPaso {p}: falta {REQUISITOS[p]}\n'
                     f'(corre antes el paso anterior)')
        print('\n' + '#' * 70 + f'\n#  PASO {p}: {titulo}\n' + '#' * 70)
        t0 = time.time()
        funcion()
        plt.close('all')
        print(f'\n[paso {p} terminado en {time.time() - t0:.0f} s]')
    print('\nListo. Resultados en analisis_galfitm/')

# -*- coding: utf-8 -*-
"""
inspeccion_visual_bd.py
=======================
Inspección visual de las descomposiciones bulbo + disco de GalfitM.

Muestra, galaxia por galaxia, las que el filtro automático acepta (por defecto las
categorías bd_ok, bulbo_no_resuelto y disco_puro de filtro_bd.py) y guarda tu decisión
con una tecla. El progreso se guarda después de CADA decisión, así que puedes cortar y
seguir otro día: al volver a abrirlo empieza por la primera galaxia sin revisar.

Figura de cada galaxia (misma escala de brillo en imagen, modelo, bulbo y disco):
  arriba : imagen | banda r suavizada (gaussiana, escala log + isofotas) | modelo B+D |
           residuo (imagen - modelo) | residuo relativo
  abajo  : bulbo  | disco | máscara de enfoque de r (suavizado sigma / suavizado 4 sigma - 1) |
           perfil radial en anillos elípticos (dato, modelo, bulbo, disco y FWHM de la PSF)
  Los paneles de r no dependen del modelo: muestran barras, brazos y anillos aunque el
  bulbo los haya absorbido (isofotas internas alargadas o giradas = posible barra).
  título : ID, categoría, tipo morfológico, B/T, n_b, Re_b, Re_d, b/a, chi2/nu

TECLAS
  a  o  1   aceptar            (visual_ok = 1)  y pasar a la siguiente
  x  o  0   rechazar           (visual_ok = 0)  y pasar a la siguiente
  d  o  2   dudosa             (visual_ok = 2)  y pasar a la siguiente
  u  o  3   una componente     (visual_ok = 3)  se describe bien con UNA sola componente
            (p. ej. una elíptica: el bulbo domina en todos los radios y el "disco" solo
            corrige el perfil). No es un mal ajuste: entra al análisis como esferoide
            con su ajuste single-Sérsic, no como bulbo + disco.
  b  barra   o  anillo/lente   v  vecino/contaminación
  s  estructura (brazos, grumos)   e  artefacto/borde   i  componentes intercambiadas
            (las etiquetas se activan/desactivan ANTES de decidir; se guardan con la decisión)
  →  /  ←   siguiente / anterior sin decidir
  n         ir a la próxima sin revisar
  f         cambiar de banda
  [  /  ]   menos / más suavizado gaussiano (paneles 'r suavizada' y 'máscara de enfoque')
  +  /  -   acercar / alejar
  q         guardar y salir

Salida: analisis_galfitm/bd/filtro/inspeccion_visual.csv
  objectId, visual (aceptada/rechazada/dudosa/una_componente), visual_ok (1/0/2/3), visual_etiquetas,
  categoria_bd, banda_vista, fecha

Para pasar la columna a las tablas del análisis:
  python3 inspeccion_visual_bd.py --unir
  agrega visual, visual_ok y visual_etiquetas a
    bd/filtro/galfitm_bd_filtrado.csv  y  bd/masa_tamano/componentes_masas.csv
  (hace una copia .bak de cada una la primera vez). La decisión vive en
  inspeccion_visual.csv, así que si vuelves a correr filtro_bd.py basta con repetir --unir.

Uso (desde analisis_conjunto/):
  python3 inspeccion_visual_bd.py                        # revisar las aceptadas
  python3 inspeccion_visual_bd.py --categorias bd_ok     # solo bd_ok
  python3 inspeccion_visual_bd.py --orden azar           # en orden aleatorio
  python3 inspeccion_visual_bd.py --revisar dudosa       # volver a mirar las dudosas
  python3 inspeccion_visual_bd.py --ids lista.csv        # solo estos IDs
  python3 inspeccion_visual_bd.py --png 30               # sin ventana: 30 PNG (para probar)
  python3 inspeccion_visual_bd.py --resumen              # cuántas llevas
  python3 inspeccion_visual_bd.py --unir
"""
import os
import sys
import glob
import argparse
import warnings
from datetime import datetime
import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib.patches import Ellipse
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import FuncFormatter, NullFormatter
from scipy.ndimage import gaussian_filter, label as nd_label
from scipy.special import gammaincinv, gamma as gamma_fn
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIG
# =============================================================================
BASE = os.path.dirname(os.path.abspath(__file__))
AG = os.path.join(BASE, 'analisis_galfitm')
FILTRO = os.path.join(AG, 'bd', 'filtro', 'galfitm_bd_filtrado.csv')
RESULT = os.path.join(AG, 'bd', 'galfitm_bd_resultados.csv')
MORFO = os.path.join(AG, 'morfologia', 'galaxias_morfologia.csv')
SS_CSV = os.path.join(AG, 'crossmatch', 'galfitm_x_lsst.csv')       # ajuste de una componente
SS_PATRONES = ['{oid}_galfitm_out.fits', '{oid}_galfitm*out*.fits', '{oid}*sersic*out*.fits']
COMP = os.path.join(AG, 'bd', 'masa_tamano', 'componentes_masas.csv')
SALIDA = os.path.join(AG, 'bd', 'filtro', 'inspeccion_visual.csv')
SALIDA_CLAUDE = os.path.join(AG, 'bd', 'filtro', 'inspeccion_visual_claude.csv')   # segundo clasificador (solo --unir --final)
PNG_DIR = os.path.join(AG, 'bd', 'filtro', 'inspeccion_png')
SUFIJO = '_galfitm_bd_out.fits'

CATEGORIAS_DEF = 'bd_ok,bulbo_no_resuelto,disco_puro'
BANDAS_PREF = ['r', 'i', 'g', 'z', 'u', 'y']
PIXSCALE = 0.2
ZP = 31.4                     # mismo punto cero que leer_output_bd.py
ZOOM_RE = 5.0                 # medio lado del recorte = ZOOM_RE x Re del disco
MIN_MEDIO_LADO = 25           # píxeles
BANDA_SUAVE = 'r'             # banda de los paneles suavizados (si falta, la que se esté viendo)
SIGMA_SUAVE = 1.5             # sigma del suavizado gaussiano [px]; se cambia con [ y ]
FACTOR_UNSHARP = 4.0          # máscara de enfoque: suavizado(sigma) / suavizado(FACTOR x sigma) - 1

DECISION = {'a': 1, '1': 1, 'x': 0, '0': 0, 'd': 2, '2': 2, 'u': 3, '3': 3}
NOMBRE = {1: 'aceptada', 0: 'rechazada', 2: 'dudosa', 3: 'una_componente'}
COLOR_DEC = {1: '#2e7d32', 0: '#c62828', 2: '#ef6c00', 3: '#6a1b9a'}
ETIQUETAS = {'b': 'barra', 'o': 'anillo', 'v': 'vecino', 's': 'estructura',
             'e': 'artefacto', 'i': 'intercambiadas'}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--categorias', default=CATEGORIAS_DEF)
ap.add_argument('--orden', choices=['id', 'azar', 'categoria'], default='categoria')
ap.add_argument('--revisar', choices=['dudosa', 'rechazada', 'aceptada', 'una_componente', 'todas'], default=None,
                help='volver a mostrar galaxias ya revisadas con esa decisión')
ap.add_argument('--ids', default=None, help='CSV o TXT con objectId: solo esas galaxias')
ap.add_argument('--png', type=int, default=0, metavar='N', help='sin ventana: guarda N figuras PNG')
ap.add_argument('--resumen', action='store_true')
ap.add_argument('--incluir-auto', action='store_true',
                help='mostrar también las que filtro_una_componente.py ya pasó a una componente')
ap.add_argument('--unir', action='store_true')
ap.add_argument('--final', action='store_true',
                help='con --unir: agrega además visual_ok_final, visual_fuente y visual_etiquetas_final '
                     '(tu decisión; si no hay, la de Claude; si no, la automática de una componente)')
ap.add_argument('--semilla', type=int, default=42)
ap.add_argument('--idioma', choices=['en', 'es'], default='en', help='textos de la figura (en = inglés)')
ap.add_argument('--limpia', action='store_true',
                help='figura para mostrar: sin contador, sin decisión ni ayuda de teclas')
args = ap.parse_args()


# =============================================================================
# TABLAS
# =============================================================================
def leer(p, **kw):
    if not os.path.exists(p):
        sys.exit(f'No encuentro {p}')
    return pd.read_csv(p, dtype={'objectId': str}, low_memory=False, **kw)


def leer_inspeccion():
    if os.path.exists(SALIDA):
        v = pd.read_csv(SALIDA, dtype={'objectId': str, 'visual_etiquetas': str})
        v['visual_etiquetas'] = v.visual_etiquetas.fillna('')
        return v.drop_duplicates('objectId', keep='last').set_index('objectId')
    return pd.DataFrame(columns=['visual', 'visual_ok', 'visual_etiquetas', 'categoria_bd',
                                 'banda_vista', 'fecha']).rename_axis('objectId')


def guardar_inspeccion(v):
    tmp = SALIDA + '.tmp'
    v.reset_index().to_csv(tmp, index=False)
    os.replace(tmp, SALIDA)           # escritura atómica: no se corrompe si se corta


COLS_DEC = ['visual', 'visual_ok', 'visual_etiquetas', 'categoria_bd', 'banda_vista', 'fecha']


def guardar_decision(oid, valores):
    """Guarda UNA decisión sin pisar el resto: vuelve a leer el CSV del disco, cambia solo la
    fila de esta galaxia y lo reescribe. Así dos ventanas abiertas a la vez (o una ventana
    abierta desde antes) no se borran las decisiones entre sí. Devuelve la tabla actualizada."""
    t = leer_inspeccion()
    t.loc[oid, COLS_DEC] = valores
    guardar_inspeccion(t)
    return t


def resumen(v, total=None):
    print(f'\nRevisadas: {len(v)}' + (f' de {total}' if total else ''))
    if len(v):
        print(v.visual.value_counts().to_string())
        et = v.visual_etiquetas[v.visual_etiquetas != ''].str.split(';').explode().value_counts()
        if len(et):
            print('etiquetas:\n' + et.to_string())
        if 'categoria_bd' in v:
            print('\n' + pd.crosstab(v.categoria_bd, v.visual).to_string())


if args.unir:
    v = leer_inspeccion()
    if v.empty:
        sys.exit('Todavía no hay inspección guardada.')
    cols = ['visual', 'visual_ok', 'visual_etiquetas']
    cols_fin = ['visual_ok_final', 'visual_fuente', 'visual_etiquetas_final']
    if args.final:
        # prioridad: humana > claude > auto (filtro_una_componente.py: una_comp_motivo no vacío -> 3)
        vc = pd.DataFrame(columns=['visual_ok', 'visual_etiquetas']).rename_axis('objectId')
        if os.path.exists(SALIDA_CLAUDE):
            vc = pd.read_csv(SALIDA_CLAUDE, dtype={'objectId': str, 'visual_etiquetas': str}) \
                .drop_duplicates('objectId', keep='last').set_index('objectId')
        else:
            print(f'(no existe {SALIDA_CLAUDE}: solo se usan tus decisiones y la clase automática)')
        fa = leer(FILTRO)
        auto = set(fa.objectId[fa.una_comp_motivo.fillna('').astype(str) != '']) if 'una_comp_motivo' in fa else set()
        disco_puro = set(fa.objectId[fa.categoria_bd == 'disco_puro']) if 'categoria_bd' in fa else set()
    for p in (FILTRO, COMP):
        if not os.path.exists(p):
            print(f'(no existe {p}, se omite)')
            continue
        if not os.path.exists(p + '.bak'):
            import shutil
            shutil.copy(p, p + '.bak')
        t = leer(p)
        t = t.drop(columns=[c for c in cols if c in t]).merge(v[cols].reset_index(), on='objectId', how='left')
        t['visual'] = t.visual.fillna('sin_revisar')
        t = t.drop(columns=[c for c in cols_fin if c in t])     # sin --final no quedan columnas finales viejas
        if args.final:
            oid = t.objectId
            en_h, en_c, en_a = oid.isin(v.index), oid.isin(vc.index), oid.isin(auto)
            t['visual_ok_final'] = np.select(
                [en_h, en_c, en_a],
                [pd.to_numeric(oid.map(v.visual_ok), errors='coerce'), pd.to_numeric(oid.map(vc.visual_ok), errors='coerce'), 3.0],
                np.nan)
            # disco_puro es una componente por definición: "aceptada" (1) ahí significa que el ajuste de
            # una componente es bueno -> 3 (solo en la columna final; visual_ok queda como se guardó)
            dp1 = oid.isin(disco_puro) & (t.visual_ok_final == 1)
            t.loc[dp1, 'visual_ok_final'] = 3.0
            t['visual_fuente'] = np.select([en_h, en_c, en_a], ['humana', 'claude', 'auto'], '')
            t['visual_etiquetas_final'] = np.where(en_h, oid.map(v.visual_etiquetas).fillna(''),
                                                   np.where(en_c, oid.map(vc.visual_etiquetas).fillna(''), ''))
        t.to_csv(p, index=False)
        print(f'{p}: {t.visual.value_counts().to_dict()}')
        if args.final:
            print(f'  disco_puro con 1 -> 3 en visual_ok_final: {int(dp1.sum())}')
            print('  final por fuente:\n' + pd.crosstab(t.visual_fuente.replace('', 'ninguna'),
                                                        t.visual_ok_final.fillna(-1).astype(int).replace(-1, 'sin')).to_string())
    sys.exit(0)

if args.resumen:
    resumen(leer_inspeccion())
    sys.exit(0)

f = leer(FILTRO)
cats = [c.strip() for c in args.categorias.split(',') if c.strip()]
# con --ids se muestran esas galaxias sea cual sea su categoría
sel = f.copy() if (args.ids or 'categoria_bd' not in f) else f[f.categoria_bd.isin(cats)].copy()
if 'una_comp_motivo' in sel and not args.incluir_auto and not args.ids:
    auto = sel.una_comp_motivo.fillna('').astype(str) != ''
    if auto.any():
        print(f'Se omiten {int(auto.sum())} galaxias que filtro_una_componente.py ya clasificó como '
              'una componente (usa --incluir-auto para verlas).')
    sel = sel[~auto]

r = leer(RESULT)
cols_r = [c for c in r.columns if c not in sel.columns or c == 'objectId']
sel = sel.merge(r[cols_r], on='objectId', how='left')
# ajuste de UNA componente (single-Sérsic): columnas SS_<param>_<banda>
if os.path.exists(SS_CSV):
    ss = leer(SS_CSV)
    params = ('MAG', 'RE', 'N', 'AR', 'PA', 'XC', 'YC')
    ren = {}
    for c in ss.columns:
        base = c.split('_')[0].upper()
        if c.upper().startswith('SS_'):
            ren[c] = c
        elif base in params or c.upper() in ('CHI2NU', 'CHI2_NU'):
            ren[c] = 'SS_' + (c.upper() if c.upper().startswith('CHI2') else c)
    ss = ss[['objectId'] + list(ren)].rename(columns=ren).drop_duplicates('objectId')
    sel = sel.merge(ss[[c for c in ss.columns if c == 'objectId' or c not in sel.columns]],
                    on='objectId', how='left')
    print('Ajuste de una componente: ' + (', '.join(sorted(c for c in ss.columns if c != 'objectId'))[:200]
          or 'no encontré columnas MAG/RE/N/AR en ' + SS_CSV))
else:
    print(f'(no encuentro {SS_CSV}: no se mostrará el ajuste de una componente)')
if 'TypeName' not in sel and os.path.exists(MORFO):
    mo = leer(MORFO)
    sel = sel.merge(mo[[c for c in ('objectId', 'TypeName', 'logM_taylor', 'z') if c in mo]],
                    on='objectId', how='left', suffixes=('', '_m'))
if args.ids:
    if args.ids.endswith('.csv'):
        ids = set(leer(args.ids).objectId)
    else:
        ids = {l.strip() for l in open(args.ids) if l.strip()}
    sel = sel[sel.objectId.isin(ids)]

v = leer_inspeccion()
if args.revisar == 'todas':
    sel = sel[sel.objectId.isin(v.index)]
elif args.revisar:
    sel = sel[sel.objectId.isin(v.index[v.visual == args.revisar])]

if args.orden == 'azar':
    sel = sel.sample(frac=1, random_state=args.semilla)
elif args.orden == 'categoria' and 'categoria_bd' in sel:
    sel = sel.assign(_o=sel.categoria_bd.map({c: i for i, c in enumerate(cats)})).sort_values(['_o', 'objectId'])
else:
    sel = sel.sort_values('objectId')
sel = sel.reset_index(drop=True)
if sel.empty:
    sys.exit('No hay galaxias para mostrar con esas opciones.')
print(f'Galaxias a mostrar: {len(sel)}  (categorías: {", ".join(cats)})')
print(f'Ya revisadas entre ellas: {int(sel.objectId.isin(v.index).sum())}')


# =============================================================================
# ARCHIVOS FITS
# =============================================================================
_indice = None


def ruta_fits(fila):
    oid = fila.objectId
    car = fila.get('carpeta') if hasattr(fila, 'get') else None
    if isinstance(car, str) and car:
        for c in (car, os.path.join(BASE, car)):
            for p in (os.path.join(c, oid + SUFIJO), os.path.join(c, oid, oid + SUFIJO)):
                if os.path.exists(p):
                    return p
    global _indice
    if _indice is None:            # búsqueda de respaldo (una sola vez)
        _indice = {}
        for patron in (os.path.join(BASE, '..', 'galaxias_*'), os.path.join(BASE, 'galaxias_*')):
            for p in glob.glob(os.path.join(patron, '**', '*' + SUFIJO), recursive=True):
                _indice[os.path.basename(p)[:-len(SUFIJO)]] = p
    return _indice.get(oid)


def ruta_fits_ss(ruta_bd, oid):
    carpeta = os.path.dirname(ruta_bd)
    for pat in SS_PATRONES:
        for p in sorted(glob.glob(os.path.join(carpeta, pat.format(oid=oid)))):
            if '_bd' not in os.path.basename(p):
                return p
    return None


def isofotas(img, xc, yc, rmax_px, n_niv=26):
    """Forma de las isofotas por momentos: para cada nivel de brillo, la región conectada al
    centro. Devuelve semieje mayor [px], elipticidad y PA (convención GALFIT: desde +y, antihorario).
    Robusto y rápido; no es un ajuste de elipses (photutils), pero muestra barras, lentes y discos."""
    h = int(min(rmax_px * 1.2, xc, yc, img.shape[1] - xc - 1, img.shape[0] - yc - 1))
    if h < 8:
        return np.empty((0, 3))
    x0, y0 = int(round(xc)), int(round(yc))
    c = np.nan_to_num(img[y0 - h:y0 + h + 1, x0 - h:x0 + h + 1])
    sm = gaussian_filter(c, 1.0)
    borde = np.r_[img[:30].ravel(), img[-30:].ravel()]
    borde = borde[np.isfinite(borde)]
    sig = 1.4826 * np.median(np.abs(borde - np.median(borde))) if len(borde) else np.nanstd(c)
    if not sig > 0:                        # borde enmascarado (todo ceros): ruido del recorte
        sig = 1.4826 * np.median(np.abs(c - np.median(c)))
    pico = sm[h, h]
    if not np.isfinite(pico) or not sig > 0 or pico <= 5 * sig:
        return np.empty((0, 3))
    yy, xx = np.indices(sm.shape)
    out = []
    for L in np.geomspace(2 * sig, 0.85 * pico, n_niv):
        lab, _ = nd_label(sm > L)
        k = lab[h, h]
        if k == 0:
            continue
        m = lab == k
        if m.sum() < 6 or m[0].any() or m[-1].any() or m[:, 0].any() or m[:, -1].any():
            continue                       # toca el borde del recorte: no es una isofota cerrada
        x, y = xx[m] - h, yy[m] - h
        cov = np.array([[np.mean(x * x), np.mean(x * y)], [np.mean(x * y), np.mean(y * y)]])
        ev, evec = np.linalg.eigh(cov)
        q = np.sqrt(max(ev[0], 1e-9) / ev[1])
        ang = np.degrees(np.arctan2(evec[1, 1], evec[0, 1]))      # eje mayor, desde +x
        pa = (ang - 90.0) % 180.0                                 # convención GALFIT
        a = np.sqrt(m.sum() / np.pi / q)
        out.append((a, 1 - q, pa))
    return np.array(sorted(out)) if out else np.empty((0, 3))


def perfil_sersic(r_px, mag, re_px, n, q):
    """Sérsic analítico (sin PSF) a lo largo del eje mayor, en cuentas/pixel."""
    n = np.clip(n, 0.2, 10)
    bn = gammaincinv(2 * n, 0.5)
    F = 10 ** (-0.4 * (mag - ZP))
    Ie = F / (2 * np.pi * q * re_px ** 2 * np.exp(bn) * n * bn ** (-2 * n) * gamma_fn(2 * n))
    return Ie * np.exp(-bn * ((r_px / re_px) ** (1 / n) - 1))


def hdu(ext, *nombres):
    up = {h.name.upper(): h for h in ext}
    for n in nombres:
        if n.upper() in up and up[n.upper()].data is not None:
            return up[n.upper()].data.astype(float)
    return None


def componente(ext, k, b):
    for h in ext:
        n = h.name.upper()
        if n.startswith(f'COMPONENT_{k}_') and n.endswith(f'_{b.upper()}') and h.data is not None:
            return h.data.astype(float)
    return None


def bandas_disponibles(ext):
    nombres = {h.name.upper() for h in ext}
    return [b for b in BANDAS_PREF if f'INPUT_{b.upper()}' in nombres]


# =============================================================================
# FIGURA
# =============================================================================
TEXTOS = {
    'en': {'img': 'image', 'mod': 'B+D model', 'res': 'residual (image − model)',
           'rel': 'residual / model', 'bul': 'bulge (Sérsic)', 'dis': 'disc (exponential)',
           'suave': '{b} band, Gaussian-smoothed ({s})',
           'unsharp': '{b} unsharp mask ({s} / {k:g}σ)',
           'nofits': 'not in the FITS file', 'p_img': 'image', 'p_mod': 'B+D model',
           'p_bul': 'bulge', 'p_dis': 'disc', 'fwhm': 'PSF FWHM',
           'p_ss': 'single Sérsic', 'sinpsf': '(no PSF)', 'una_c': 'Single Sérsic',
           'eps': 'ellipticity', 'iso_x': 'isophote semi-major axis [arcsec]', 'sin_iso': 'no closed isophotes',
           'iso_t': 'isophote shape (r): lines = 1 − b/a of bulge (red), disc (blue), single Sérsic (green)',
           'elipse_ss': 'green dotted ellipse: single-Sérsic R$_e$',
           'xlab': 'semi-major axis of the disc ellipse [arcsec]',
           'ejes': 'offset from centre [arcsec]', 'tipo': 'type', 'banda': 'band',
           'elipses': 'dashed ellipses: R$_e$ of bulge (red) and disc (blue)',
           'contador': '{i}/{n}, inspected {h}', 'dec': 'saved decision', 'etq': 'flags',
           'dec_n': {1: 'accepted', 0: 'rejected', 2: 'uncertain', 3: 'single component'},
           'ayuda': ('a/1 accept · x/0 reject · d/2 uncertain · u/3 single component · b o v s e i flags · '
                     '←/→ move · n next uninspected · f band · +/- zoom · [/] smoothing · q quit')},
    'es': {'img': 'imagen', 'mod': 'modelo B+D', 'res': 'residuo (imagen − modelo)',
           'rel': 'residuo / modelo', 'bul': 'bulbo (Sérsic)', 'dis': 'disco (exponencial)',
           'suave': 'banda {b} suavizada ({s})',
           'unsharp': 'máscara de enfoque {b} ({s} / {k:g}σ)',
           'nofits': 'no está en el FITS', 'p_img': 'imagen', 'p_mod': 'modelo B+D',
           'p_bul': 'bulbo', 'p_dis': 'disco', 'fwhm': 'FWHM PSF',
           'p_ss': 'una componente', 'sinpsf': '(sin PSF)', 'una_c': 'Una componente',
           'eps': 'elipticidad', 'iso_x': 'semieje mayor de la isofota [arcsec]', 'sin_iso': 'sin isofotas cerradas',
           'iso_t': 'forma de las isofotas (r): líneas = 1 − b/a de bulbo (rojo), disco (azul), una comp. (verde)',
           'elipse_ss': 'elipse verde punteada: R$_e$ de una componente',
           'xlab': 'semieje mayor de la elipse del disco [arcsec]',
           'ejes': 'distancia al centro [arcsec]', 'tipo': 'tipo', 'banda': 'banda',
           'elipses': 'elipses: R$_e$ del bulbo (rojo) y del disco (azul)',
           'contador': '{i}/{n}, revisadas {h}', 'dec': 'decisión guardada', 'etq': 'etiquetas',
           'dec_n': NOMBRE,
           'ayuda': ('a/1 aceptar · x/0 rechazar · d/2 dudosa · u/3 una componente · b o v s e i etiquetas · '
                     '←/→ mover · n próxima sin revisar · f banda · +/- zoom · [/] suavizado · q salir')},
}
T = TEXTOS[args.idioma]
ETQ_EN = {'barra': 'bar', 'anillo': 'ring/lens', 'vecino': 'neighbour', 'estructura': 'structure',
          'artefacto': 'artefact', 'intercambiadas': 'swapped components'}
def num(fila, *cols):
    for c in cols:
        if c in fila and pd.notna(fila[c]):
            try:
                return float(fila[c])
            except (TypeError, ValueError):
                pass
    return np.nan


def norm_asinh(img, lims, stretch=0.03):
    v1, v2 = lims
    x = np.clip((img - v1) / (v2 - v1 + 1e-30), 0, 1)
    return np.arcsinh(x / stretch) / np.arcsinh(1 / stretch)


def perfil(img, xc, yc, q, pa_deg, rmax):
    """Mediana en anillos elípticos (PA de GALFIT: desde +y, antihorario)."""
    yy, xx = np.indices(img.shape)
    dx, dy = xx - xc, yy - yc
    t = np.deg2rad(pa_deg + 90.0)
    xr = dx * np.cos(t) + dy * np.sin(t)
    yr = -dx * np.sin(t) + dy * np.cos(t)
    r = np.hypot(xr, yr / max(q, 0.05))
    bordes = np.unique(np.r_[0, np.geomspace(1, max(rmax, 2), 22)])
    rc, val = [], []
    for a, b_ in zip(bordes[:-1], bordes[1:]):
        m = (r >= a) & (r < b_) & np.isfinite(img)
        if m.sum() >= 3:
            rc.append(0.5 * (a + b_))
            val.append(np.median(img[m]))
    return np.array(rc), np.array(val)


def mu(I):
    I = np.asarray(I, float)
    out = np.full(I.shape, np.nan)
    ok = I > 0
    out[ok] = -2.5 * np.log10(I[ok] / PIXSCALE ** 2) + ZP
    return out


def dibujar(fig, fila, banda_pref=None, zoom=1.0, estado=None, sigma=SIGMA_SUAVE):
    fig.clf()
    p = ruta_fits(fila)
    if p is None:
        fig.text(0.5, 0.5, f'{fila.objectId}: no encuentro {fila.objectId}{SUFIJO}', ha='center')
        return None, []
    with fits.open(p) as ext:
        bandas = bandas_disponibles(ext)
        if not bandas:
            fig.text(0.5, 0.5, f'{fila.objectId}: el FITS no tiene INPUT_<banda>', ha='center')
            return None, []
        b = banda_pref if banda_pref in bandas else bandas[0]
        B = b.upper()
        img = hdu(ext, f'INPUT_{B}')
        mod = hdu(ext, f'MODEL_{B}')
        res = hdu(ext, f'RESIDUAL_{B}')
        bul = componente(ext, 1, b)
        dis = componente(ext, 2, b)
        b_s = BANDA_SUAVE if BANDA_SUAVE in bandas else b
        img_s = hdu(ext, f'INPUT_{b_s.upper()}')
    if res is None and img is not None and mod is not None:
        res = img - mod
    mod_ss = None
    p_ss = ruta_fits_ss(p, fila.objectId)
    if p_ss:
        try:
            with fits.open(p_ss) as e2:
                mod_ss = hdu(e2, f'MODEL_{B}')
        except Exception:
            mod_ss = None

    xc = num(fila, f'D_XC_{b}', f'B_XC_{b}', 'D_XC_r', 'B_XC_r')
    yc = num(fila, f'D_YC_{b}', f'B_YC_{b}', 'D_YC_r', 'B_YC_r')
    if not np.isfinite(xc) or not np.isfinite(yc):
        yc, xc = (np.array(img.shape) - 1) / 2 + 1
    xc, yc = xc - 1, yc - 1                              # GALFIT cuenta desde 1
    re_d = num(fila, f'D_RE_{b}', 'D_RE_r', 'D_RE_i')
    if not np.isfinite(re_d):
        re_d = 1.678 * num(fila, f'D_RS_{b}', 'D_RS_r')
    re_b = num(fila, f'B_RE_{b}', 'B_RE_r', 'B_RE_i')
    q_d = num(fila, f'D_AR_{b}', 'D_AR_r')
    pa_d = num(fila, f'D_PA_{b}', 'D_PA_r')
    q_b = num(fila, f'B_AR_{b}', 'B_AR_r')
    pa_b = num(fila, f'B_PA_{b}', 'B_PA_r')
    medio = max(MIN_MEDIO_LADO, ZOOM_RE * (re_d if np.isfinite(re_d) else 10)) / zoom
    x0, x1 = int(max(0, xc - medio)), int(min(img.shape[1], xc + medio + 1))
    y0, y1 = int(max(0, yc - medio)), int(min(img.shape[0], yc + medio + 1))
    cut = lambda a: None if a is None else a[y0:y1, x0:x1]

    # escala común: del ruido del fondo al pico del modelo (stretch asinh fuerte: se ve el disco externo)
    sig_bg = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res))) if res is not None else 0
    lims = (-1.0 * sig_bg, np.nanpercentile(cut(mod), 99.9))
    ext_ = [x0 - xc, x1 - xc, y0 - yc, y1 - yc]
    ext_ = [e * PIXSCALE for e in ext_]                  # ejes en arcsec respecto del centro
    # suavizado gaussiano y máscara de enfoque (pasa-banda) de la banda r
    lim_s = lims
    if img_s is not None:
        f_img = np.where(np.isfinite(img_s), img_s, 0.0)
        suave = gaussian_filter(f_img, sigma)
        ancho = gaussian_filter(f_img, FACTOR_UNSHARP * sigma)
        sig_s = 1.4826 * np.nanmedian(np.abs(suave - np.nanmedian(suave)))
        unsharp = np.where(ancho > 3 * max(sig_s, 1e-30), suave / ancho - 1, np.nan)
        lim_s = (-1.0 * sig_s, np.nanmax(cut(suave)))
    else:
        suave = unsharp = None
    gs = fig.add_gridspec(2, 5, height_ratios=[1, 1])
    s_txt = f'σ = {sigma * PIXSCALE:.2f}"'
    paneles = [(cut(img), T['img'], 'gray_r', 'lin'),
               (cut(suave), T['suave'].format(b=b_s, s=s_txt), 'gray_r', 'suave'),
               (cut(mod), T['mod'], 'gray_r', 'lin'), (cut(res), T['res'], 'RdBu_r', 'res'),
               (None, T['rel'], 'RdBu_r', 'rel'),
               (cut(bul), T['bul'], 'Reds', 'lin'), (cut(dis), T['dis'], 'Blues', 'lin'),
               (cut(unsharp), T['unsharp'].format(b=b_s, s=s_txt, k=FACTOR_UNSHARP), 'PuOr_r', 'um')]
    pos = [(0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (1, 0), (1, 1), (1, 2)]
    sig = np.nanstd(cut(res)) if res is not None else 1
    for (d, tit, cmap, modo), (i, j) in zip(paneles, pos):
        ax = fig.add_subplot(gs[i, j])
        if modo == 'rel':
            if res is not None and mod is not None:
                rel = cut(res) / np.where(cut(mod) > 3 * sig_bg, cut(mod), np.nan)
                im = ax.imshow(rel, cmap=cmap, vmin=-0.5, vmax=0.5, origin='lower', extent=ext_)
                cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, ticks=[-0.5, -0.25, 0, 0.25, 0.5])
                cb.ax.tick_params(labelsize=7)
        elif d is None:
            ax.text(0.5, 0.5, T['nofits'], ha='center', va='center', transform=ax.transAxes)
        elif modo == 'suave':
            # escala logarítmica + isofotas: una barra o lente se ve como isofotas internas
            # alargadas o giradas respecto de las externas
            lo_, hi_ = 2 * max(sig_s, 1e-30), np.nanmax(d)
            dl = np.log10(np.clip(d, lo_, None))
            ax.imshow(dl, cmap=cmap, vmin=np.log10(lo_), vmax=np.log10(hi_), origin='lower', extent=ext_)
            if hi_ > 4 * lo_:
                niv = np.geomspace(3 * sig_s, 0.8 * hi_, 12)
                ax.contour(d, levels=niv, colors='#ff7f0e', linewidths=0.6, origin='lower', extent=ext_)
        elif modo == 'um':
            # escala fijada fuera del núcleo (el pico central domina la máscara de enfoque)
            yy_, xx_ = np.indices(d.shape)
            rr_ = np.hypot(xx_ - (xc - x0), yy_ - (yc - y0))
            fuera = np.isfinite(d) & (rr_ > 1.5 * FACTOR_UNSHARP * sigma)
            v = np.nanpercentile(np.abs(d[fuera]), 90) if fuera.sum() > 20 else 0.1
            vmax_ = max(np.nanmax(np.abs(d)), 2 * v)
            # simétrica-logarítmica: lineal hasta v (brazos, barra) y comprimida en el núcleo
            im = ax.imshow(d, cmap=cmap, origin='lower', extent=ext_,
                           norm=SymLogNorm(linthresh=v, linscale=1.5, vmin=-vmax_, vmax=vmax_))
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.ax.tick_params(labelsize=7)
        elif modo == 'res':
            im = ax.imshow(d, cmap=cmap, vmin=-4 * sig, vmax=4 * sig, origin='lower', extent=ext_)
            cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02, ticks=[-4 * sig, 0, 4 * sig])
            cb.ax.set_yticklabels(['−4σ', '0', '+4σ'], fontsize=7)
        else:
            ax.imshow(norm_asinh(d, lims), cmap=cmap, vmin=0, vmax=1, origin='lower', extent=ext_)
        for re_, q__, pa__, col in ((re_b, q_b, pa_b, '#c62828'), (re_d, q_d, pa_d, '#1565c0')):
            if np.isfinite(re_) and modo != 'rel':
                ax.add_patch(Ellipse((0, 0), 2 * re_ * PIXSCALE,
                                     2 * re_ * PIXSCALE * (q__ if np.isfinite(q__) else 1),
                                     angle=(pa__ if np.isfinite(pa__) else 0) + 90,
                                     fill=False, color=col, lw=1.0, ls='--'))
        if tit == T['img']:
            re_s = num(fila, f'SS_RE_{b}', 'SS_RE', 'SS_RE_r')
            q_s = num(fila, f'SS_AR_{b}', 'SS_AR', 'SS_AR_r')
            pa_s = num(fila, f'SS_PA_{b}', 'SS_PA', 'SS_PA_r')
            if np.isfinite(re_s):
                ax.add_patch(Ellipse((0, 0), 2 * re_s * PIXSCALE, 2 * re_s * PIXSCALE * (q_s if np.isfinite(q_s) else 1),
                                     angle=(pa_s if np.isfinite(pa_s) else 0) + 90,
                                     fill=False, color='#2e7d32', lw=1.2, ls=':'))
        ax.set_title(tit, fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(False)
        if j == 0:
            ax.set_ylabel(T['ejes'], fontsize=7)

    ax = fig.add_subplot(gs[1, 3])
    rmax = medio
    q_ = q_d if np.isfinite(q_d) else 1.0
    pa_ = pa_d if np.isfinite(pa_d) else 0.0
    mu_dato = []
    for d, lab, col, ls in ((img, T['p_img'], 'k', ''), (mod, T['p_mod'], '0.4', '-'),
                            (bul, T['p_bul'], '#c62828', '--'), (dis, T['p_dis'], '#1565c0', '--')):
        if d is None:
            continue
        rr, vv = perfil(d, xc, yc, q_, pa_, rmax)
        if d is img:
            ax.plot(rr * PIXSCALE, mu(vv), 'o', color=col, ms=3.5, label=lab)
            mu_dato = mu(vv)
        else:
            ax.plot(rr * PIXSCALE, mu(vv), ls, color=col, lw=1.6, label=lab)
    # ajuste de una componente: modelo del FITS (con PSF) o Sérsic analítico (sin PSF)
    if mod_ss is not None:
        rr, vv = perfil(mod_ss, xc, yc, q_, pa_, rmax)
        ax.plot(rr * PIXSCALE, mu(vv), ':', color='#2e7d32', lw=2.2, label=T['p_ss'])
    else:
        ss_n_, ss_re_, ss_q_, ss_mag_ = (num(fila, f'SS_N_{b}', 'SS_N', 'SS_N_r'), num(fila, f'SS_RE_{b}', 'SS_RE', 'SS_RE_r'),
                                         num(fila, f'SS_AR_{b}', 'SS_AR', 'SS_AR_r'), num(fila, f'SS_MAG_{b}', 'SS_MAG', 'SS_MAG_r'))
        if all(np.isfinite([ss_n_, ss_re_, ss_mag_])):
            rr = np.geomspace(0.5, rmax, 80)
            q_s = ss_q_ if np.isfinite(ss_q_) else 1.0
            # el perfil de los datos se mide en anillos con la forma del disco: se compara en el eje mayor
            ax.plot(rr * PIXSCALE, mu(perfil_sersic(rr, ss_mag_, ss_re_, ss_n_, q_s)), ':', color='#2e7d32',
                    lw=2.2, label=T['p_ss'] + ' ' + T['sinpsf'])
    fw = num(fila, f'FWHM_{b}', 'FWHM_PIX', 'fwhm_pix', 'FWHM_REF')
    if np.isfinite(fw):
        ax.axvline(fw * PIXSCALE, color='0.5', ls=':', lw=1, label=T['fwhm'])
    mu_ok = np.asarray(mu_dato)[np.isfinite(mu_dato)]
    if len(mu_ok):                          # eje y en el rango de los datos (no lo estira una componente)
        ax.set_ylim(mu_ok.min() - 1.0, mu_ok.max() + 2.0)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.set_xlabel(T['xlab'], fontsize=8)
    ax.set_ylabel(r'$\mu$ [mag/arcsec$^2$]', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6.5, loc='upper right')

    # ---------- forma de las isofotas (banda r) ----------
    ax = fig.add_subplot(gs[1, 4])
    base_iso = img_s if img_s is not None else img
    iso = isofotas(base_iso, xc, yc, medio) if base_iso is not None else np.empty((0, 3))
    ax2 = ax.twinx()
    if len(iso) and np.isfinite(fw):
        iso = iso[iso[:, 0] >= 0.5 * fw]        # dentro de media FWHM la forma no está definida
    if len(iso):
        a_iso = iso[:, 0] * PIXSCALE
        ref = pa_d if np.isfinite(pa_d) else np.median(iso[:, 2])
        pa_rel = (iso[:, 2] - ref + 90) % 180 - 90 + ref          # PA cerca del del disco (sin saltos de 180)
        ax.plot(a_iso, iso[:, 1], 'o-', color='k', ms=3.5, lw=1.2, label=T['eps'])
        ax2.plot(a_iso, pa_rel, 's--', color='#8e24aa', ms=3, lw=1, label='PA')
        ss_q_iso = num(fila, f'SS_AR_{b}', 'SS_AR', 'SS_AR_r')
        for qq, pp, col in ((q_b, pa_b, '#c62828'), (q_d, pa_d, '#1565c0'), (ss_q_iso, num(fila, f'SS_PA_{b}', 'SS_PA', 'SS_PA_r'), '#2e7d32')):
            if np.isfinite(qq):
                ax.axhline(1 - qq, color=col, lw=1.2, ls='-', alpha=0.8)
            if np.isfinite(pp):
                ax2.axhline((pp - ref + 90) % 180 - 90 + ref, color=col, lw=0.8, ls=':', alpha=0.8)
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(0, max(0.7, np.nanmax(iso[:, 1]) + 0.05))
        lo_pa, hi_pa = np.nanmin(pa_rel), np.nanmax(pa_rel)
        ax2.set_ylim(min(lo_pa, ref - 30) - 5, max(hi_pa, ref + 30) + 5)
    else:
        ax.text(0.5, 0.5, T['sin_iso'], ha='center', va='center', transform=ax.transAxes, fontsize=8)
    if np.isfinite(fw):
        ax.axvline(fw * PIXSCALE, color='0.5', ls=':', lw=1)
    ax.set_xlabel(T['iso_x'], fontsize=8)
    ax.set_ylabel(T['eps'], fontsize=8)
    ax2.set_ylabel('PA [deg]', fontsize=8, color='#8e24aa')
    ax.tick_params(labelsize=7)
    ax2.tick_params(labelsize=7, colors='#8e24aa')
    ax.set_title(T['iso_t'], fontsize=7.5)

    def ss(par):
        return num(fila, f'SS_{par}_{b}', f'SS_{par}', f'SS_{par}_r', f'SS_{par}_i')
    ss_n, ss_re, ss_q, ss_mag = ss('N'), ss('RE'), ss('AR'), ss('MAG')
    ss_chi = num(fila, 'SS_CHI2NU', 'SS_chi2nu', f'SS_CHI2NU_{b}')
    bt = num(fila, f'BT_{b}', 'BT', 'BT_r')
    n_b = num(fila, f'B_N_{b}', 'B_N', 'B_N_r')
    chi = num(fila, 'chi2nu', 'CHI2NU')
    tipo = fila.get('TypeName', '?') if hasattr(fila, 'get') else '?'
    linea1 = f'{fila.objectId}   {T["tipo"]}: {tipo}   {T["banda"]} {b}'
    if not args.limpia:
        linea1 += (f'   [{fila.get("categoria_bd", "?")}]   ('
                   + T['contador'].format(i=estado['i'] + 1, n=estado['n'], h=estado['hechas']) + ')')
    linea2 = (f'B/T = {bt:.2f}   $n_b$ = {n_b:.2f}   $R_{{e,b}}$ = {re_b * PIXSCALE:.2f}"'
              f'   $R_{{e,d}}$ = {re_d * PIXSCALE:.2f}"   $(b/a)_b$ = {q_b:.2f}   $(b/a)_d$ = {q_d:.2f}'
              f'   $\\chi^2_\\nu$ = {chi:.2f}      ' + T['elipses'])
    linea3 = ''
    if np.isfinite(ss_re):
        dchi = f'   Δ$\\chi^2_\\nu$ (1c − B+D) = {ss_chi - chi:+.3f}' if np.isfinite(ss_chi) and np.isfinite(chi) else ''
        linea3 = ('\n' + T['una_c'] + f':  n = {ss_n:.2f}   $R_e$ = {ss_re * PIXSCALE:.2f}"   b/a = {ss_q:.2f}'
                  + (f'   $\\chi^2_\\nu$ = {ss_chi:.2f}' if np.isfinite(ss_chi) else '') + dchi
                  + '      ' + T['elipse_ss'])
    fig.suptitle(linea1 + '\n' + linea2 + linea3, fontsize=10, y=0.998)
    if not args.limpia:
        previa = estado.get('previa')
        etq = estado.get('etiquetas', set())
        etq_txt = [ETQ_EN.get(e, e) if args.idioma == 'en' else e for e in sorted(etq)]
        txt = f'{T["dec"]}: ' + (T['dec_n'][previa] if previa is not None else '—')
        txt += f'     {T["etq"]}: ' + (', '.join(etq_txt) if etq_txt else '—')
        fig.text(0.01, 0.005, txt, fontsize=10, color=COLOR_DEC.get(previa, 'k'), fontweight='bold')
        fig.text(0.99, 0.005, T['ayuda'], fontsize=8, ha='right', color='0.35')
    fig.subplots_adjust(left=0.04, right=0.97, top=0.87, bottom=0.07, wspace=0.32, hspace=0.3)
    return b, bandas


# =============================================================================
# MODO PNG (sin ventana)
# =============================================================================
if args.png:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    os.makedirs(PNG_DIR, exist_ok=True)
    fig = plt.figure(figsize=(20, 8.5))
    for i, fila in sel.head(args.png).iterrows():
        prev = v.loc[fila.objectId] if fila.objectId in v.index else None
        est = {'i': i, 'n': len(sel), 'hechas': len(v),
               'previa': None if prev is None else int(prev.visual_ok),
               'etiquetas': set() if prev is None else set(filter(None, str(prev.visual_etiquetas).split(';')))}
        dibujar(fig, fila, estado=est)
        fig.savefig(os.path.join(PNG_DIR, f'{fila.objectId}_insp.png'), dpi=90)
    print(f'{min(args.png, len(sel))} figuras en {PNG_DIR}/')
    sys.exit(0)

# =============================================================================
# MODO INTERACTIVO
# =============================================================================
import matplotlib
import matplotlib.pyplot as plt
NO_INTERACTIVOS = ('agg', 'pdf', 'svg', 'ps', 'cairo', 'template')


def es_interactivo():
    bk = matplotlib.get_backend().lower()
    return not (bk in NO_INTERACTIVOS or 'inline' in bk)


if not es_interactivo():
    errores = []
    for bk in ('TkAgg', 'QtAgg', 'Qt5Agg', 'GTK3Agg'):
        try:
            plt.switch_backend(bk)
            plt.figure().clf(); plt.close('all')      # fuerza a abrir la conexión gráfica
            break
        except Exception as e:
            errores.append(f'{bk}: {type(e).__name__}: {str(e)[:80]}')
            plt.switch_backend('Agg')
    if not es_interactivo():
        print('\nNo pude abrir una ventana: matplotlib no encontró un backend gráfico.')
        print('\n'.join('  ' + e for e in errores))
        print('\nSoluciones (cualquiera):')
        print('  sudo apt install python3-tk          (si usas el python del sistema)')
        print('  conda install tk                      (si usas conda / (base))')
        print('  pip install pyqt5')
        print('Mientras tanto puedes mirar figuras sin ventana:  python3 inspeccion_visual_bd.py --png 30')
        sys.exit(1)
print(f'Ventana con backend {matplotlib.get_backend()}. Haz clic en la ventana antes de usar las teclas.')
for k in list(plt.rcParams):
    if k.startswith('keymap.'):
        plt.rcParams[k] = []              # que matplotlib no use las teclas (s, o, f, q, ...)

st = {'i': 0, 'banda': None, 'zoom': 1.0, 'etiquetas': set(), 'sigma': SIGMA_SUAVE}
pendientes = ~sel.objectId.isin(v.index)
if args.revisar is None and pendientes.any():
    st['i'] = int(np.argmax(pendientes.to_numpy()))
fig = plt.figure(figsize=(20, 8.5))


def cargar_previa():
    oid = sel.objectId.iloc[st['i']]
    if oid in v.index:
        r_ = v.loc[oid]
        st['previa'] = int(r_.visual_ok)
        st['etiquetas'] = set(filter(None, str(r_.visual_etiquetas).split(';')))
    else:
        st['previa'] = None
        st['etiquetas'] = set()


def refrescar():
    est = {'i': st['i'], 'n': len(sel), 'hechas': int(sel.objectId.isin(v.index).sum()),
           'previa': st['previa'], 'etiquetas': st['etiquetas']}
    b, bandas = dibujar(fig, sel.iloc[st['i']], st['banda'], st['zoom'], est, st['sigma'])
    st['banda_vista'], st['bandas'] = b, bandas
    fig.canvas.draw_idle()


def mover(paso):
    st['i'] = int(np.clip(st['i'] + paso, 0, len(sel) - 1))
    st['zoom'] = 1.0
    cargar_previa()
    refrescar()


def on_key(ev):
    global v
    k = (ev.key or '').lower()
    oid = sel.objectId.iloc[st['i']]
    if k in DECISION:
        val = DECISION[k]
        v = guardar_decision(oid, [
            NOMBRE[val], val, ';'.join(sorted(st['etiquetas'])),
            sel.categoria_bd.iloc[st['i']] if 'categoria_bd' in sel else '',
            st.get('banda_vista', ''), datetime.now().isoformat(timespec='seconds')])
        print(f'{oid}: {NOMBRE[val]} {sorted(st["etiquetas"]) or ""}')
        if st['i'] == len(sel) - 1:
            print('Era la última galaxia.')
            refrescar()
            return
        mover(1)
    elif k in ETIQUETAS:
        st['etiquetas'] ^= {ETIQUETAS[k]}
        refrescar()
    elif k in ('right', 'down', ' '):
        mover(1)
    elif k in ('left', 'up', 'backspace'):
        mover(-1)
    elif k == 'n':
        resto = (~sel.objectId.isin(v.index)).to_numpy()
        prox = np.where(resto[st['i'] + 1:])[0]
        if len(prox):
            mover(int(prox[0]) + 1)
        else:
            print('No quedan galaxias sin revisar después de esta.')
    elif k == 'f':
        bs = st.get('bandas') or []
        if bs:
            j = (bs.index(st['banda_vista']) + 1) % len(bs) if st.get('banda_vista') in bs else 0
            st['banda'] = bs[j]
            refrescar()
    elif k == '[':
        st['sigma'] = max(0.5, st['sigma'] / 1.5)
        refrescar()
    elif k == ']':
        st['sigma'] = min(10.0, st['sigma'] * 1.5)
        refrescar()
    elif k in ('+', '='):
        st['zoom'] *= 1.5
        refrescar()
    elif k in ('-', '_'):
        st['zoom'] /= 1.5
        refrescar()
    elif k == 'q':
        plt.close(fig)


fig.canvas.mpl_connect('key_press_event', on_key)
cargar_previa()
refrescar()
plt.show()
v = leer_inspeccion()                  # cada decisión ya se guardó; aquí solo se relee (no se pisa el disco)
resumen(v[v.index.isin(sel.objectId)], total=len(sel))
print(f'\nGuardado en {SALIDA}')
print('Para pasarlo a las tablas:  python3 inspeccion_visual_bd.py --unir')

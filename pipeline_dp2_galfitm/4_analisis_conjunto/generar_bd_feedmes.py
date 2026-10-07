# -*- coding: utf-8 -*-
"""
generar_bd_feedmes.py
=====================
Genera feedmes Bulge + Disc para GalfitM a partir del ajuste single-Sersic,
para la muestra confiable A + B (procesar_muestra.py). Adaptado de
make_bd_feedmes_V2.py (DP1).

Receta de valores iniciales: Vika et al. (2014, MNRAS 444, 3603), sec. 2.2.2
(INICIO = 'vika2014'); alternativas 'v2' (tu V2) y 'proporcional'.
Receta original de V2:
  Bulbo (sersic):  m = m_ss + 0.75,  Re = 0.5 * Re_ss,  n_init = 4
  Disco (expdisk): m = m_ss + 0.65,  Rs = Re_ss / 1.678
  Ambos comparten centro, q y PA iniciales (los del single-Sersic).

Cambios respecto a V2 (necesarios para DP2):
  1. BANDAS POR GALAXIA. En DP2 cada galaxia tiene su propio juego de bandas
     (p. ej. g,r,i,z,y o i,z,y). Las bandas se leen de la linea A1) de su
     feedme single-Sersic, en vez de exigir ugrizy (que dejaba fuera a casi todas).
  2. BANDAS MALAS FUERA DEL AJUSTE. Las bandas que no pasaron los cortes de
     calidad (imagen vacia, fondo malo, flag de GalfitM...) se quitan de las
     lineas A, A1, A2, C, D, F y J. Como GalfitM ajusta todas las bandas juntas,
     una banda mala sesgaria tambien a las demas.
  3. FLAGS POLINOMICOS <= NUMERO DE BANDAS. Con 4 bandas, "cuadratico" (3)
     casi equivale a libre por banda; se avisa y se limita a nb.
  4. SE CORRE DENTRO DE CADA CARPETA. Las rutas del feedme quedan relativas a la
     carpeta de la galaxia, asi los archivos galfit.NN no se mezclan y se puede
     correr en paralelo (ver correr_bd.sh).
  5. LA SECCION DE COMPONENTES SE REESCRIBE COMPLETA (bulbo + disco + cielo), en
     vez de buscar un bloque '#---Galaxy' con una expresion regular.

Entrada:
  analisis_galfitm/muestra/muestra_AB.csv     (procesar_muestra.py)
  galfitm_inputs/<ID>/<ID>.feedme              (feedme single-Sersic = template)
Salida:
  galfitm_inputs/<ID>/<ID>_bd.feedme
  lista_bd.txt                                 (carpetas a correr con correr_bd.sh)
  resumen_bd_feedmes.csv                       (bandas usadas y excluidas por galaxia)

Uso:
    python3 generar_bd_feedmes.py
    ./correr_bd.sh 4          # 4 galaxias en paralelo
"""

import os
import re
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# =============================================================================
# CONFIG
# =============================================================================
BASE        = os.path.dirname(os.path.abspath(__file__))
MUESTRA     = os.path.join(BASE, 'analisis_galfitm', 'muestra', 'muestra_AB.csv')
GALFITM_DIR = os.path.join(BASE, 'galfitm_inputs')
SUFIJO      = '_bd'
MIN_BANDAS  = 3              # minimo de bandas buenas para intentar el ajuste
BANDAS_OBLIGATORIAS = ['g', 'r', 'i']   # control extra: sin estas tres no se genera el feedme
NO_SOBRESCRIBIR = True   # si <ID>_bd.feedme ya existe se deja igual (protege ajustes
                         # en curso o ya hechos; pon False para regenerarlos todos)

EXCLUIR_BANDAS_MALAS = True  # False = usar todas las bandas del template

# --- Valores iniciales -------------------------------------------------------
# 'vika2014' (por defecto): receta EXACTA de Vika et al. (2014, MNRAS 444, 3603),
#     seccion 2.2.2, para descomposicion multibanda con GalfitM:
#       bulbo : m = m_ss + 0.75, Re = 0.5 Re_ss, n = n_ss, b/a = 0.8, PA = 10 deg
#       disco : m = m_ss + 0.65, Re = Re_ss (Rs = Re/1.678), n = 1 fijo,
#               b/a y PA del single-Sersic
#       Re, n, b/a, PA y centro CONSTANTES con lambda; magnitudes libres por banda;
#       mismo centro para bulbo y disco; cielo FIJO.
#       Si un parametro varia con lambda en el single-Sersic, se usa la mediana.
#     Restricciones (archivo G): mag 5-35, Re 0.04-600", n 0.1-15, centro comun
#     y desplazamiento <= sqrt(s/8) px respecto al single-Sersic.
# 'v2': tu make_bd_feedmes_V2.py (misma receta de radios y flujos que Vika 2014,
#     pero n_b = 4, b/a y PA del bulbo = disco, Re/n cuadraticos, sin restricciones).
# 'proporcional': Re_b = 0.25 Re_ss y B/T segun n_ss (alternativa propia, sin paper).
INICIO = 'vika2014'

# Flags polinomicos por receta (1 = constante en lambda, 2 = lineal,
# 3 = cuadratico, nb = libre por banda). Se limitan a nb.
FLAGS_POR_RECETA = {
    'vika2014':     dict(bulge_re=1, bulge_n=1, bulge_ar=1, bulge_pa=1,
                         disc_rs=1, disc_ar=1, disc_pa=1),
    'v2':           dict(bulge_re=3, bulge_n=3, bulge_ar=1, bulge_pa=1,
                         disc_rs=3, disc_ar=1, disc_pa=1),
    'proporcional': dict(bulge_re=1, bulge_n=1, bulge_ar=1, bulge_pa=1,
                         disc_rs=1, disc_ar=1, disc_pa=1),
}
FLAGS = FLAGS_POR_RECETA[INICIO]
CIELO_LIBRE = (INICIO == 'v2')   # Vika 2014: cielo fijo
USAR_RESTRICCIONES = (INICIO != 'v2')

DISC_COMPONENT = 'expdisk'   # 'expdisk' o 'sersic' (n=1 fijo)
PIXSCALE = 0.2               # arcsec/pixel (para las restricciones de Re)

# Receta proporcional (alternativa)
DISC_RE_FRAC    = 1.1
BULGE_RE_FRAC_P = 0.25
PSF_FWHM_PIX    = 3.5
BULGE_Q_MIN     = 0.6
BT_POR_N = [(1.5, 0.15, 2.0), (3.0, 0.30, 3.0), (np.inf, 0.50, 4.0)]

# Lineas del header con un valor por banda (se recortan si se excluyen bandas)
LINEAS_POR_BANDA = ('A', 'A1', 'A2', 'C', 'D', 'F', 'J')


# =============================================================================
# FUNCIONES
# =============================================================================
def leer_template(path):
    """Separa el header (lineas A..W) de los componentes (desde el primer ' 0) ')."""
    lineas = Path(path).read_text().splitlines(keepends=True)
    fin = len(lineas)
    for k, l in enumerate(lineas):
        if re.match(r'^\s*0\)\s', l):
            fin = k
            break
    # cortar despues de la ultima linea de control (A..W): se descartan los
    # comentarios del componente viejo ('# Component number: 1', etc.)
    claves = [k for k in range(fin) if re.match(r'^\s*[A-Z]\d?\)\s', lineas[k])]
    return lineas[:claves[-1] + 1] if claves else lineas[:fin]


def valor_linea(header, clave):
    """Devuelve el primer token de la linea 'clave)'."""
    for l in header:
        m = re.match(rf'^\s*{re.escape(clave)}\)\s+(\S+)', l)
        if m:
            return m.group(1)
    return None


def recortar_header(header, n_total, indices):
    """Deja solo las bandas `indices` en las lineas con un valor por banda."""
    out = []
    for l in header:
        m = re.match(r'^(\s*)([A-Z]\d?)\)(\s+)(\S+)(.*)$', l.rstrip('\n'))
        if m and m.group(2) in LINEAS_POR_BANDA:
            vals = m.group(4).split(',')
            if len(vals) == n_total:
                vals = [vals[i] for i in indices]
                l = f'{m.group(1)}{m.group(2)}){m.group(3)}{",".join(vals)}{m.group(5)}\n'
        out.append(l)
    return out


def ajustar_header(header, oid, out_fits):
    """Salida B), imagenes de componentes en W) y rutas relativas a la carpeta."""
    txt = ''.join(header)
    txt = re.sub(r'^(\s*B\)\s+)\S+', rf'\g<1>{out_fits}', txt, flags=re.MULTILINE)
    if re.search(r'^\s*W\)', txt, flags=re.MULTILINE):
        txt = re.sub(r'^(\s*W\)\s+)\S+', r'\g<1>input,model,residual,component',
                     txt, flags=re.MULTILINE)
    else:
        txt += 'W) input,model,residual,component\n'
    # 'galfitm_inputs/<ID>/archivo.fits' -> 'archivo.fits' (se corre dentro de la carpeta)
    txt = re.sub(rf'[^\s,]*{oid}/', '', txt)
    return txt


def make_bd_components(row, bands, flags, disc_component=DISC_COMPONENT):
    """Bloques bulbo + disco + cielo segun la receta INICIO."""
    nb = len(bands)
    g = lambda p: np.array([float(row[f'{p}_{b}']) for b in bands])
    mag, re_, ar, pa, xc, yc, n_ = g('MAG'), g('RE'), g('AR'), g('PA'), g('XC'), g('YC'), g('N')
    sky = [float(row.get(f'SKY_{b}')) if pd.notna(row.get(f'SKY_{b}')) else 0.0 for b in bands]
    x0, y0 = np.median(xc), np.median(yc)
    fl = {k: min(v, nb) for k, v in flags.items()}
    cte = lambda v: [float(np.median(v))] * nb      # parametro constante con lambda

    if INICIO == 'vika2014':
        mag_b, mag_d = mag + 0.75, mag + 0.65
        re_b, re_d = cte(0.5 * re_), cte(re_)
        n_b = cte(n_)
        ar_b, pa_b = [0.8] * nb, [10.0] * nb
        ar_d, pa_d = cte(ar), cte(pa)
    elif INICIO == 'v2':
        mag_b, mag_d = mag + 0.75, mag + 0.65
        re_b, re_d = list(0.5 * re_), list(re_)
        n_b = [4.0] * nb
        ar_b, pa_b, ar_d, pa_d = list(ar), list(pa), list(ar), list(pa)
    else:  # proporcional
        n_ss = float(np.median(n_))
        bt, n0 = next((bt, nn) for lim, bt, nn in BT_POR_N if n_ss < lim)
        mag_b, mag_d = mag - 2.5 * np.log10(bt), mag - 2.5 * np.log10(1 - bt)
        re_d = cte(DISC_RE_FRAC * re_)
        re_b = [max(BULGE_RE_FRAC_P * float(np.median(re_)), PSF_FWHM_PIX)] * nb
        n_b = [n0] * nb
        ar_b, pa_b = [max(float(np.median(ar)), BULGE_Q_MIN)] * nb, cte(pa)
        ar_d, pa_d = cte(ar), cte(pa)
    size_d = [r / 1.678 for r in re_d] if disc_component == 'expdisk' else list(re_d)

    fmt = lambda vals, p=4: ','.join(f'{v:.{p}f}' for v in vals)
    bulge = (f"# ---------- Bulge (Sersic)  [receta: {INICIO}] ----------\n"
             " 0) sersic\n"
             f" 1) {fmt([x0]*nb)} 1\n"
             f" 2) {fmt([y0]*nb)} 1\n"
             f" 3) {fmt(mag_b)} {nb}\n"
             f" 4) {fmt(re_b)} {fl['bulge_re']}\n"
             f" 5) {fmt(n_b)} {fl['bulge_n']}\n"
             f" 9) {fmt(ar_b)} {fl['bulge_ar']}\n"
             f"10) {fmt(pa_b)} {fl['bulge_pa']}\n"
             " Z) 0\n")
    if disc_component == 'expdisk':
        disc = ("# ---------- Disc (expdisk) ----------\n"
                " 0) expdisk\n"
                f" 1) {fmt([x0]*nb)} 1\n"
                f" 2) {fmt([y0]*nb)} 1\n"
                f" 3) {fmt(mag_d)} {nb}\n"
                f" 4) {fmt(size_d)} {fl['disc_rs']}\n"
                f" 9) {fmt(ar_d)} {fl['disc_ar']}\n"
                f"10) {fmt(pa_d)} {fl['disc_pa']}\n"
                " Z) 0\n")
    else:
        disc = ("# ---------- Disc (Sersic, n=1 fijo) ----------\n"
                " 0) sersic\n"
                f" 1) {fmt([x0]*nb)} 1\n"
                f" 2) {fmt([y0]*nb)} 1\n"
                f" 3) {fmt(mag_d)} {nb}\n"
                f" 4) {fmt(size_d)} {fl['disc_rs']}\n"
                f" 5) {fmt([1.0]*nb)} 0\n"
                f" 9) {fmt(ar_d)} {fl['disc_ar']}\n"
                f"10) {fmt(pa_d)} {fl['disc_pa']}\n"
                " Z) 0\n")
    cielo = ("# ---------- Sky ----------\n"
             " 0) sky\n"
             f" 1) {fmt(sky, 3)} {nb if CIELO_LIBRE else 0}\n"
             f" 2) {fmt([0.0]*nb, 1)} 0\n"
             f" 3) {fmt([0.0]*nb, 1)} 0\n"
             " Z) 0\n")
    return bulge + "\n" + disc + "\n" + cielo


def restricciones_vika(n_pix_lado):
    """Archivo de restricciones de Vika et al. (2014), seccion 2.2.2.
    Componente 1 = bulbo, 2 = disco. Sintaxis de GALFIT:
      'a to b'  rango absoluto;  'dmin dmax'  rango relativo al valor inicial;
      '1_2 x offset'  mantiene el mismo centro en las dos componentes."""
    dxy = np.sqrt(n_pix_lado / 8.0)
    re_min, re_max = 0.04 / PIXSCALE, 600.0 / PIXSCALE
    return ("# Restricciones de Vika et al. (2014, MNRAS 444, 3603), sec. 2.2.2\n"
            "# componente  parametro  restriccion\n"
            "  1_2         x          offset\n"
            "  1_2         y          offset\n"
            f"  1           x          {-dxy:.2f} {dxy:.2f}\n"
            f"  1           y          {-dxy:.2f} {dxy:.2f}\n"
            "  1           mag        5 to 35\n"
            "  2           mag        5 to 35\n"
            f"  1           re         {re_min:.2f} to {re_max:.1f}\n"
            f"  2           rs         {re_min / 1.678:.2f} to {re_max / 1.678:.1f}\n"
            "  1           n          0.1 to 15\n")


# =============================================================================
# MAIN
# =============================================================================
if __name__ == '__main__':
    if not os.path.exists(MUESTRA):
        sys.exit(f'No encuentro {MUESTRA}. Corre antes procesar_muestra.py')
    df = pd.read_csv(MUESTRA, dtype={'objectId': str}, low_memory=False)
    print(f'Muestra A + B: {len(df)} galaxias')
    if 'categoria' in df:
        print('  ' + str(df.categoria.value_counts().to_dict()))

    avisos_flag = set()
    filas, lista, errores = [], [], []
    for _, row in df.iterrows():
        oid = row.objectId
        ruta = row.get('ruta_carpeta')
        carpeta = os.path.join(ruta if isinstance(ruta, str) and ruta else BASE,
                               'galfitm_inputs', oid)
        template = os.path.join(carpeta, f'{oid}.feedme')
        if not os.path.exists(template):
            errores.append((oid, 'sin feedme single-Sersic'))
            continue
        if NO_SOBRESCRIBIR and os.path.exists(os.path.join(carpeta, f'{oid}{SUFIJO}.feedme')):
            lista.append(carpeta)
            filas.append({'objectId': oid, 'bandas_template': '', 'bandas_bd': '',
                          'excluidas': '', 'n_bandas': np.nan, 'ya_existia': True})
            continue
        try:
            header = leer_template(template)
            etiquetas = valor_linea(header, 'A1')
            if etiquetas is None:
                raise ValueError('el template no tiene linea A1) con las bandas')
            bandas_tpl = [b.strip() for b in etiquetas.split(',')]

            # bandas buenas = las que tienen todos los parametros tras los cortes
            def buena(b):
                return all(pd.notna(row.get(f'{p}_{b}'))
                           for p in ('MAG', 'RE', 'AR', 'PA', 'XC', 'YC'))
            usar = [i for i, b in enumerate(bandas_tpl) if buena(b) or not EXCLUIR_BANDAS_MALAS]
            fuera = [b for i, b in enumerate(bandas_tpl) if i not in usar]
            bandas = [bandas_tpl[i] for i in usar]
            if not EXCLUIR_BANDAS_MALAS:
                # valores iniciales de una banda mala: mediana de las buenas
                for b in bandas:
                    if not buena(b):
                        for p in ('MAG', 'RE', 'AR', 'PA', 'XC', 'YC'):
                            vals = [row[f'{p}_{x}'] for x in bandas_tpl if buena(x)]
                            row[f'{p}_{b}'] = np.median(vals)
            faltan = [b for b in BANDAS_OBLIGATORIAS if b not in bandas]
            if faltan:
                raise ValueError(f'faltan bandas obligatorias: {",".join(faltan)}')
            if len(bandas) < MIN_BANDAS:
                raise ValueError(f'solo {len(bandas)} bandas buenas ({",".join(bandas)})')

            for k, v in FLAGS.items():
                if v > 1 and v >= len(bandas):
                    avisos_flag.add((k, v, len(bandas)))

            header = recortar_header(header, len(bandas_tpl), usar)
            out_fits = f'{oid}_galfitm{SUFIJO}_out.fits'
            cabecera = ajustar_header(header, oid, out_fits)
            if USAR_RESTRICCIONES:
                m = re.search(r'^\s*H\)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', cabecera, re.MULTILINE)
                lado = (int(m.group(2)) - int(m.group(1)) + 1) if m else 601
                fcons = f'{oid}{SUFIJO}.constraints'
                Path(carpeta, fcons).write_text(restricciones_vika(lado))
                cabecera = re.sub(r'^(\s*G\)\s+)\S+', rf'\g<1>{fcons}', cabecera, flags=re.MULTILINE)
            texto = cabecera + '\n' + make_bd_components(row, bandas, FLAGS)
            Path(carpeta, f'{oid}{SUFIJO}.feedme').write_text(texto)

            lista.append(carpeta)
            filas.append({'objectId': oid, 'bandas_template': ','.join(bandas_tpl),
                          'bandas_bd': ','.join(bandas), 'excluidas': ','.join(fuera),
                          'n_bandas': len(bandas)})
        except Exception as e:
            errores.append((oid, str(e)))

    Path(BASE, 'lista_bd.txt').write_text('\n'.join(lista) + '\n')
    res = pd.DataFrame(filas)
    res.to_csv(os.path.join(BASE, 'resumen_bd_feedmes.csv'), index=False)

    print(f'\nFeedmes bulbo + disco generados: {len(lista)}')
    print(f'Errores / sin template          : {len(errores)}')
    for oid, msg in errores[:5]:
        print(f'   {oid}: {msg}')
    if len(res):
        if 'ya_existia' in res:
            print(f'Feedmes que ya existian (no se tocaron): {int(res.ya_existia.fillna(False).sum())}')
        print(f'Numero de bandas por galaxia    : {res.n_bandas.dropna().astype(int).value_counts().sort_index().to_dict()}')
        con_fuera = (res.excluidas != '').sum()
        print(f'Galaxias con bandas excluidas   : {con_fuera}')
    for k, v, nb in sorted(avisos_flag):
        print(f'AVISO: {k}={v} con galaxias de {nb} bandas -> {min(v, nb)} coeficientes = '
              f'libre por banda (sin suavidad en lambda)')
    print(f'\nLista para correr: lista_bd.txt\nSiguiente paso:    ./correr_bd.sh <N procesos>')

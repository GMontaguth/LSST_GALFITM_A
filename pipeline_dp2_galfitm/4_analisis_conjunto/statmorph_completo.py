#!/usr/bin/env python3
"""
Morfología con statmorph-lsst: TODOS los parámetros (no paramétricos +
Sérsic simple + Sérsic doble), con la PSF de cada banda.

Lee las mismas carpetas que usas para GalfitM:
    BASE/<objectId>/<objectId>_sci_<b>.fits
    BASE/<objectId>/<objectId>_sig_<b>.fits
    BASE/<objectId>/<objectId>_psf_<b>.fits   (para los ajustes Sérsic)
    BASE/<objectId>/<objectId>_mask.fits      (opcional, máscara SEP: 1 = enmascarado)

y escribe un CSV con una fila por (galaxia, banda), solo en g, r, i.
La columna gri_ok marca las galaxias con las tres bandas medidas y flag <= 1.

Uso:
    python3 statmorph_completo.py /ruta/galaxias_todas            # todas
    python3 statmorph_completo.py /ruta/galaxias_todas --ncpu 8   # en paralelo
    python3 statmorph_completo.py /ruta/galaxias_todas --png 20   # + figuras QA de las 20 primeras
    python3 statmorph_completo.py /ruta/galaxias_todas --lista muestra.csv   # solo los objectId de la tabla
"""
import os
import sys
import glob
import argparse
import warnings
import traceback
from multiprocessing import Pool

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel, convolve
from photutils.segmentation import detect_sources
from scipy import ndimage as ndi
from astropy.utils import lazyproperty

# --- Encontrar statmorph-lsst aunque no este instalado con pip -------------
# Se busca en: variable STATMORPH_LSST_DIR, carpeta statmorph-lsst junto a este
# script, o una o dos carpetas mas arriba (p. ej. Galfitm/statmorph-lsst).
try:
    import statmorph_lsst
except ImportError:
    _aqui = os.path.dirname(os.path.abspath(__file__))
    _candidatos = [os.environ.get('STATMORPH_LSST_DIR', ''),
                   os.path.join(_aqui, 'statmorph-lsst'),
                   os.path.join(_aqui, '..', 'statmorph-lsst'),
                   os.path.join(_aqui, '..', '..', 'statmorph-lsst')]
    for _c in _candidatos:
        for _r in (_c, os.path.join(_c, 'src')):       # repo con o sin carpeta src/
            if _c and os.path.isdir(os.path.join(_r, 'statmorph_lsst')):
                sys.path.insert(0, os.path.abspath(_r))
                break
        else:
            continue
        break
    try:
        import statmorph_lsst
    except ImportError:
        sys.exit('No encuentro statmorph_lsst. Instalalo con\n'
                 '   pip install -e "/ruta/a/statmorph-lsst"\n'
                 'o define STATMORPH_LSST_DIR con la ruta de la carpeta statmorph-lsst.')
from statmorph_lsst import SourceMorphology, _quantity_names

# ======================= CONFIGURACIÓN =======================
BANDAS        = ['g', 'r', 'i']   # las mismas obligatorias de procesar_muestra.py
PIXEL_SCALE   = 0.2          # arcsec/pix (coadds LSST)
ZP_NJY        = 31.4         # AB: m = -2.5 log10(f[nJy]) + 31.4
SB_ISOFOTAS   = [24.0, 25.0, 26.0]   # mag/arcsec² para la asimetría isofotal
NSIGMA_SEG    = 1.5          # umbral de detección para el segmap (en sigma)
NPIX_SEG      = 10           # área mínima del segmento
SMOOTH_FWHM   = 3.0          # px, suavizado antes de detectar
CUTOUT_HALF   = 150          # px; recorte alrededor del centro (None = stamp completo)
NITER_BH_MID  = 100          # basin-hopping del multimode (100 recomendado en producción)
PSF_HALF      = 25           # px; la PSF se recorta a 51x51 para acelerar la convolución
# =============================================================

# Parche: con pandas >= 3 la interpolación de la máscara devuelve un arreglo de
# solo lectura y statmorph-lsst falla ('assignment destination is read-only').
import statmorph_lsst.statmorph_lsst as _sml
_interp_orig = _sml._interpolate_missing_pixels
def _interp_copia(*a, **k):
    img, m = _interp_orig(*a, **k)
    return np.array(img, copy=True), m
_sml._interpolate_missing_pixels = _interp_copia

DOBLE_Q = list(_sml._doublesersic_quantity_names)

# Radios en píxeles que convertimos también a arcsec
RADIOS = ['rpetro_circ', 'rpetro_ellip', 'rmax_circ', 'rmax_ellip',
          'rhalf_circ', 'rhalf_ellip', 'r20', 'r50', 'r80',
          'sersic_rhalf', 'doublesersic_rhalf1', 'doublesersic_rhalf2']


class SourceMorphologyCorr(SourceMorphology):
    """statmorph-lsst completo, con la corrección de la smoothness."""

    @lazyproperty
    def _sky_smoothness(self):
        # Corrección: statmorph-lsst 0.7.1 normaliza por sum(cielo), que es ~0
        # en una imagen sin cielo y hace explotar S. Lotz+04 y el statmorph
        # original usan el residuo positivo POR PÍXEL, que es lo que la
        # fórmula S = (ap_diff - area*sky)/ap_flux necesita.
        bkg = self._cutout_stamp_maskzeroed[self._slice_skybox]
        if bkg.size == 0:
            return -99.0
        bkg_diff = self._sky_smoothness_residual.copy()
        bkg_diff[~np.isfinite(bkg_diff) | (bkg_diff < 0)] = 0.0
        return np.sum(bkg_diff) / float(bkg.size)



def leer_psf(fpsf):
    """PSF recortada, sin valores negativos/NaN y normalizada."""
    psf = fits.getdata(fpsf).astype(np.float64)
    psf[~np.isfinite(psf) | (psf < 0)] = 0.0
    ny, nx = psf.shape
    cy, cx = np.unravel_index(np.argmax(psf), psf.shape)
    h = min(PSF_HALF, cy, cx, ny - 1 - cy, nx - 1 - cx)
    psf = psf[cy - h:cy + h + 1, cx - h:cx + h + 1]
    return psf / psf.sum()


def sb_a_flujo_pix(sb):
    """mag/arcsec² -> nJy/pix."""
    return 10 ** ((ZP_NJY - sb) / 2.5) * PIXEL_SCALE ** 2


def recortar(arr, cy, cx, half):
    if half is None:
        return arr, 0, 0
    ny, nx = arr.shape
    y0, y1 = max(0, cy - half), min(ny, cy + half + 1)
    x0, x1 = max(0, cx - half), min(nx, cx + half + 1)
    return arr[y0:y1, x0:x1], y0, x0


def preparar(sci, sig, mask_sep):
    """Devuelve imagen sin cielo, sigma, máscara booleana y segmap (galaxia = 1)."""
    ny, nx = sci.shape
    cy, cx = ny // 2, nx // 2

    bad = ~np.isfinite(sci) | ~np.isfinite(sig) | (sig <= 0)
    if mask_sep is not None:
        bad |= mask_sep.astype(bool)

    sci, y0, x0 = recortar(sci, cy, cx, CUTOUT_HALF)
    sig, _, _ = recortar(sig, cy, cx, CUTOUT_HALF)
    bad, _, _ = recortar(bad, cy, cx, CUTOUT_HALF)
    cy, cx = cy - y0, cx - x0

    img = np.where(bad, 0.0, sci).astype(np.float64)
    sig = np.where(bad, np.nanmedian(sig[~bad]), sig).astype(np.float64)

    # 1) Detección preliminar para enmascarar fuentes y medir el cielo
    kern = Gaussian2DKernel(SMOOTH_FWHM / 2.355)
    _, med0, std0 = sigma_clipped_stats(img[~bad], sigma=3.0)
    sm = convolve(img - med0, kern, mask=bad)
    seg0 = detect_sources(sm, NSIGMA_SEG * std0, NPIX_SEG)
    fuentes = np.zeros_like(bad) if seg0 is None else \
        ndi.binary_dilation(seg0.data > 0, iterations=5)
    _, sky, _ = sigma_clipped_stats(img[~bad & ~fuentes], sigma=3.0)

    img = img - sky
    img[bad] = 0.0

    # 2) Segmap final sobre la imagen sin cielo, umbral por píxel con el sigma
    sm = convolve(img, kern, mask=bad)
    seg = detect_sources(sm, NSIGMA_SEG * sig, NPIX_SEG)
    if seg is None:
        raise RuntimeError('no se detectó ninguna fuente')
    lab = seg.data[cy, cx]
    if lab == 0:  # centro no detectado: segmento más cercano
        yy, xx = np.nonzero(seg.data)
        k = np.argmin((yy - cy) ** 2 + (xx - cx) ** 2)
        lab = seg.data[yy[k], xx[k]]
    segmap = (seg.data == lab).astype(np.int32)

    # Otros segmentos (vecinos que SEP no tapó) -> máscara
    otros = (seg.data > 0) & (seg.data != lab)
    mask = bad | ndi.binary_dilation(otros, iterations=2) & (segmap == 0)

    return img, sig, mask, segmap, sky


def medir_galaxia(tarea):
    warnings.simplefilter("ignore")
    gdir, oid, qa_dir, doble = tarea
    filas = []
    mf = os.path.join(gdir, f'{oid}_mask.fits')
    mask_sep = fits.getdata(mf) if os.path.exists(mf) else None
    iso_levels = [sb_a_flujo_pix(s) for s in SB_ISOFOTAS]

    for b in BANDAS:
        fsci = os.path.join(gdir, f'{oid}_sci_{b}.fits')
        fsig = os.path.join(gdir, f'{oid}_sig_{b}.fits')
        fila = {'objectId': oid, 'band': b}
        if not (os.path.exists(fsci) and os.path.exists(fsig)):
            fila['flag'] = 5            # banda obligatoria sin imagen
            fila['error'] = 'falta_banda'
            filas.append(fila)
            continue
        try:
            sci = fits.getdata(fsci).astype(np.float64)
            sig = fits.getdata(fsig).astype(np.float64)
            img, sig, mask, segmap, sky = preparar(sci, sig, mask_sep)
            fila['sky_subtracted'] = sky
            fila['mask_frac'] = mask.mean()
            fpsf = os.path.join(gdir, f'{oid}_psf_{b}.fits')
            psf = leer_psf(fpsf) if os.path.exists(fpsf) else None
            fila['psf_usada'] = psf is not None

            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                m = SourceMorphologyCorr(
                    img, segmap, 1, mask=mask, weightmap=sig, psf=psf,
                    asymmetry_isophotes=iso_levels,
                    niter_bh_mid=NITER_BH_MID,
                    include_doublesersic=doble)

            for q in list(_quantity_names) + (DOBLE_Q if doble else []):
                v = getattr(m, q)
                if q == 'isophote_asymmetry':
                    vals = v if isinstance(v, (list, tuple, np.ndarray)) else [v] * len(SB_ISOFOTAS)
                    for s, vv in zip(SB_ISOFOTAS, vals):
                        fila[f'a_iso_{s:g}'] = vv
                else:
                    fila[q] = v
            for r in RADIOS:
                v = fila.get(r, -99.0)
                fila[f'{r}_arcsec'] = v * PIXEL_SCALE if v > 0 else -99.0
            fila['flag'] = m.flag
            fila['flag_sersic'] = m.flag_sersic
            fila['sersic_runtime'] = getattr(m, 'sersic_runtime', -99.0)
            if doble:
                fila['flag_doublesersic'] = m.flag_doublesersic
            # ángulos en grados (statmorph los da en radianes, antihorario desde +x)
            for q in ['sersic_theta', 'doublesersic_theta1', 'doublesersic_theta2']:
                if q in fila and fila[q] != -99.0:
                    fila[q + '_deg'] = np.degrees(fila[q])
            fila['runtime'] = m.runtime
            if qa_dir:
                fila['_qa'] = (img, mask, segmap, m)
                figura_qa(fila, os.path.join(qa_dir, f'{oid}_{b}.png'))
                fila.pop('_qa')
        except Exception as e:
            fila['flag'] = 4
            fila['error'] = repr(e)[:200]
        filas.append(fila)
    return filas


def figura_qa(fila, outpng):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from statmorph_lsst.utils.image_diagnostics import make_figure
    m = fila['_qa'][3]
    try:
        fig = make_figure(m)
        fig.savefig(outpng, dpi=80); plt.close(fig)
        return
    except Exception:
        pass   # si falla, la figura simple de abajo
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    img, mask, segmap, m = fila['_qa']
    v1, v2 = np.percentile(img[~mask], [1, 99.7])
    norm = lambda x: np.arcsinh(np.clip((x - v1) / (v2 - v1), 0, 1) / 0.05)
    xa, ya = m._asymmetry_center
    xa += m.xmin_stamp; ya += m.ymin_stamp

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.5))
    ax[0].imshow(norm(img), cmap='gray', origin='lower')
    ax[0].add_patch(Circle((xa, ya), m.rpetro_circ, fill=False, color='cyan', label='r_petro'))
    ax[0].add_patch(Circle((xa, ya), m.r50, fill=False, color='yellow', label='r50'))
    ax[0].plot(xa, ya, '+r'); ax[0].legend(loc='upper right', fontsize=8)
    ax[0].set_title(f"{fila['objectId']}  {fila['band']}")
    ax[1].imshow(norm(img), cmap='gray', origin='lower')
    ax[1].imshow(np.ma.masked_where(~mask, mask), cmap='Reds', alpha=0.5, origin='lower')
    ax[1].contour(segmap, levels=[0.5], colors='lime', linewidths=1)
    ax[1].set_title('mask (red), segmap (green)')
    ax[2].axis('off')
    txt = (f"flag = {m.flag}\nC = {m.concentration:.3f}\nA = {m.asymmetry:.3f}\n"
           f"S = {m.smoothness:.3f}\nG = {m.gini:.3f}\nM20 = {m.m20:.3f}\n"
           f"M = {m.multimode:.3f}  I = {m.intensity:.3f}  D = {m.deviation:.3f}\n"
           f"A_shape = {m.shape_asymmetry:.3f}\nA_outer = {m.outer_asymmetry:.3f}\n"
           f"substructure = {m.substructure:.3f}\n"
           f"r_petro = {m.rpetro_circ*PIXEL_SCALE:.2f}\"  r50 = {m.r50*PIXEL_SCALE:.2f}\"\n"
           f"<S/N>/pix = {m.sn_per_pixel:.2f}")
    ax[2].text(0.02, 0.98, txt, va='top', family='monospace', fontsize=10)
    for a in ax[:2]: a.axis('off')
    plt.tight_layout(); plt.savefig(outpng, dpi=110); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('base')
    ap.add_argument('--out', default='statmorph_completo.csv')
    ap.add_argument('--ncpu', type=int, default=1)
    ap.add_argument('--png', type=int, default=0, help='figuras QA para las N primeras galaxias')
    ap.add_argument('--lista', default=None,
                    help='tabla (CSV o parquet) con los objectId a medir; si no se da, mide todo')
    ap.add_argument('--col', default='objectId', help='columna de la tabla con el ID')
    ap.add_argument('--sin-doble', action='store_true',
                    help='no ajustar el Sérsic doble (mucho más rápido)')
    args = ap.parse_args()

    if not os.path.isdir(args.base):
        sys.exit(f'No existe la carpeta: {args.base}')

    # Busca <ID>_sci_<banda>.fits a cualquier profundidad (sirve si las
    # galaxias están directamente en base/ o dentro de subcarpetas de lotes)
    galaxias = {}
    for raiz, _, archivos in os.walk(args.base):
        for f in archivos:
            if '_sci_' in f and f.endswith('.fits'):
                oid = f.split('_sci_')[0]
                if oid.isdigit():
                    galaxias.setdefault((raiz, oid), True)
    gdirs = sorted(galaxias)
    print(f'{len(gdirs)} galaxias encontradas en {args.base}')

    if args.lista:
        if not os.path.isfile(args.lista):
            sys.exit(f'No encuentro la lista: {os.path.abspath(args.lista)}\n'
                     f'Revisa que el archivo exista o pasa la ruta completa entre comillas.')
        if args.lista.endswith('.parquet'):
            tab = pd.read_parquet(args.lista, columns=[args.col])
        else:
            tab = pd.read_csv(args.lista, usecols=[args.col], dtype={args.col: str})
        ids = set(tab[args.col].dropna().astype(str).str.strip().str.replace(r'\.0$', '', regex=True))
        en_disco = {oid for _, oid in gdirs}
        gdirs = [g for g in gdirs if g[1] in ids]
        faltan = sorted(ids - en_disco)
        print(f'Lista {os.path.basename(args.lista)}: {len(ids)} IDs -> '
              f'{len(gdirs)} se miden, {len(faltan)} no están en disco')
        if faltan:
            fn = os.path.splitext(args.out)[0] + '_no_encontradas.txt'
            with open(fn, 'w') as fh:
                fh.write('\n'.join(faltan) + '\n')
            print(f'   IDs no encontrados guardados en {fn}')

    if not gdirs:
        print('\nNo encontré archivos <ID>_sci_<banda>.fits. Contenido de la carpeta:')
        for x in sorted(os.listdir(args.base))[:15]:
            print('   ', x)
        sys.exit(1)

    filas = []
    qa_dir = os.path.join(os.path.dirname(os.path.abspath(args.out)), 'qa_statmorph')
    if args.png: os.makedirs(qa_dir, exist_ok=True)

    with Pool(args.ncpu) as pool:
        tareas = [(d, oid, qa_dir if k < args.png else None, not args.sin_doble)
                  for k, (d, oid) in enumerate(gdirs)]
        for i, res in enumerate(pool.imap(medir_galaxia, tareas), 1):
            filas.extend(res)
            # Si las primeras galaxias fallan TODAS, hay un problema de
            # instalación o de datos: parar en vez de seguir horas sin medir nada
            if i == min(10, len(tareas)):
                errs = [f.get('error') for f in filas if f.get('flag') == 4]
                medidas = [f for f in filas if f.get('flag', 4) <= 3]
                if not medidas and errs:
                    pool.terminate()
                    sys.exit('\nFallaron todas las mediciones de las primeras galaxias. '
                             f'Error:\n   {errs[0]}')
            if i % 10 == 0 or i == len(gdirs):
                print(f'  {i}/{len(gdirs)}')
                pd.DataFrame(filas).to_csv(args.out, index=False)   # guardado parcial

    df = pd.DataFrame(filas)
    if df.empty or 'flag' not in df:
        sys.exit('No se midió ninguna galaxia.')
    # Galaxia completa = las tres bandas medidas con flag <= 1
    ok = df.assign(_ok=df['flag'] <= 1).groupby('objectId')['_ok'].transform('all')
    df['gri_ok'] = ok
    df.to_csv(args.out, index=False)
    print(f'\nListo: {args.out}  ({len(df)} filas)')
    print('flag por banda (5 = falta la banda):')
    print(df.groupby('band')['flag'].value_counts().unstack(fill_value=0))
    n_ok = df.loc[df.gri_ok, 'objectId'].nunique()
    print(f"\nGalaxias con g, r e i buenas (flag <= 1): {n_ok} de {df.objectId.nunique()}")


if __name__ == '__main__':
    main()

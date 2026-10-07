# -*- coding: utf-8 -*-
"""
actualizar_todo.py
==================
Un solo comando para mantener todo al dia mientras se van sumando carpetas.
Recorre las carpetas galaxias_* y corre SOLO lo que falta:

  1. En cada carpeta: ajuste single-Sersic de las galaxias sin salida y
     leer_output_galfitm.py si hay salidas nuevas.
  2. En analisis_conjunto/: une las tablas y rehace los cortes sobre TODAS las
     galaxias (las metricas de imagen se guardan; solo se miden las nuevas).
     Con --con-clasificacion ademas corrige por extincion y K y clasifica.
  3. Bulbo + disco y statmorph SOLO para las galaxias que pasan los cortes
     (categorias A y B) y que todavia no tienen resultado.
  4. Junta los resultados de statmorph en una tabla.
  5. Avisa de galaxias que ya tienen bulbo+disco o statmorph pero que con los
     cortes actuales quedaron fuera (y opcionalmente borra esos archivos).

Uso (desde analisis_conjunto/, con los scripts, calc_kcor_fun.py, TU
statmorph_completo.py y el ejecutable de GalfitM en esta carpeta).
seleccionar_discos.py NO se usa.

IMPORTANTE: si ya hay un correr_bd.sh corriendo en otra terminal, usa --sin-bd
hasta que termine, para no ajustar dos veces la misma galaxia a la vez.

    python3 actualizar_todo.py --listar        # solo muestra que falta, no corre nada
    python3 actualizar_todo.py                 # corre todo lo que falta (3 procesos)
    python3 actualizar_todo.py --nproc 4
    python3 actualizar_todo.py --sin-ss        # no lanza ajustes single-Sersic
    python3 actualizar_todo.py --sin-bd --sin-statmorph
    python3 actualizar_todo.py --solo-statmorph      # solo statmorph de las A+B que faltan
    python3 actualizar_todo.py --limpiar-sobrantes   # borra resultados de galaxias que salieron
"""

import os
import sys
import glob
import shutil
import argparse
import subprocess
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
EXE = 'galfitm-1.4.4-linux-x86_64'
CATEGORIAS_TRABAJO = ['A', 'B']        # a quienes se les corre bulbo+disco y statmorph
MUESTRA_AB = os.path.join(BASE, 'analisis_galfitm', 'muestra', 'muestra_clasificada.csv')
SALIDA_SM  = os.path.join(BASE, 'analisis_galfitm', 'statmorph_resultados.csv')
DIR_SM     = os.path.join(BASE, 'analisis_galfitm', 'statmorph_corridas')

# statmorph: se usa TU script statmorph_completo.py (statmorph-lsst), en esta carpeta.
# Se le pasa, por carpeta, una lista con SOLO las galaxias que faltan.
STATMORPH_SCRIPT = 'statmorph_completo.py'
STATMORPH_ARGS   = ['--sin-doble', '--png', '0']   # sin Sersic doble y sin figuras

ap = argparse.ArgumentParser()
ap.add_argument('carpetas', nargs='*', help='carpetas galaxias_* (por defecto ../galaxias_*)')
ap.add_argument('--nproc', type=int, default=3)
ap.add_argument('--listar', action='store_true', help='solo mostrar lo pendiente')
ap.add_argument('--sin-ss', action='store_true', help='no correr single-Sersic')
ap.add_argument('--sin-bd', action='store_true', help='no correr bulbo+disco')
ap.add_argument('--sin-statmorph', action='store_true', help='no correr statmorph')
ap.add_argument('--solo-statmorph', action='store_true',
                help='solo statmorph de las A+B que faltan, usando los cortes ya calculados')
ap.add_argument('--reintentar-fallidas', action='store_true',
                help='volver a correr el single-Sersic de las galaxias que ya fallaron')
ap.add_argument('--con-clasificacion', action='store_true',
                help='ademas correr corregir_ext_kcor.py y clasificar_morfologia.py')
ap.add_argument('--limpiar-sobrantes', action='store_true',
                help='borrar bulbo+disco/statmorph de galaxias fuera de A+B')
args = ap.parse_args()
SALTAR_PREVIOS = args.solo_statmorph      # no tocar single-Sersic ni rehacer los cortes
if args.solo_statmorph:
    args.sin_bd = True

# busca galaxias_* al lado de esta carpeta (analisis_conjunto/../) o dentro de ella
carpetas = args.carpetas or sorted(set(
    glob.glob(os.path.join(BASE, '..', 'galaxias_*')) + glob.glob(os.path.join(BASE, 'galaxias_*'))))
carpetas = sorted(set(os.path.abspath(c) for c in carpetas
                      if os.path.isdir(os.path.join(c, 'galfitm_inputs'))))
if not carpetas:
    sys.exit(f'No encontre carpetas galaxias_*/galfitm_inputs ni en {BASE} ni en la carpeta de arriba')


NECESARIOS = ['combinar_carpetas.py', 'procesar_muestra.py']
if args.solo_statmorph:
    NECESARIOS = ['statmorph_completo.py']
if not args.listar:
    faltan = [f for f in NECESARIOS if not os.path.exists(os.path.join(BASE, f))]
    if faltan:
        sys.exit(f'Faltan en {BASE}: {", ".join(faltan)}\n'
                 'Deben estar en la misma carpeta que actualizar_todo.py.')


def correr(cmd, cwd, titulo):
    print(f'\n>>> {titulo}', flush=True)
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        sys.exit(f'Fallo: {" ".join(cmd)} (en {cwd})')


def copiar_si_falta(nombre, destino):
    src, dst = os.path.join(BASE, nombre), os.path.join(destino, nombre)
    if not os.path.exists(dst) and os.path.exists(src):
        shutil.copy2(src, dst)
        if nombre.endswith('.sh') or nombre == EXE:
            os.chmod(dst, 0o755)


def ruta_galaxia(fila):
    return os.path.join(fila.ruta_carpeta, 'galfitm_inputs', fila.objectId)


def tiene_bd(d, oid):
    return os.path.exists(os.path.join(d, f'{oid}_galfitm_bd_out.fits'))


def ids_statmorph_hechos():
    """objectId que ya aparecen en alguna corrida de statmorph."""
    hechos = set()
    for f in glob.glob(os.path.join(DIR_SM, '*.csv')) + [SALIDA_SM]:
        if os.path.basename(f).startswith('pendientes_'):
            continue                                  # listas de entrada, no resultados
        if os.path.exists(f):
            try:
                hechos |= set(pd.read_csv(f, usecols=['objectId'], dtype={'objectId': str}).objectId)
            except Exception:
                pass
    return hechos


SM_HECHOS = ids_statmorph_hechos()


def tiene_sm(d, oid):
    return oid in SM_HECHOS


# =============================================================================
# 1. SINGLE-SERSIC EN CADA CARPETA
# =============================================================================
print('=' * 72)
print(f'{"carpeta":22s} {"galaxias":>9s} {"con ajuste":>11s} {"pendientes":>11s} '
      f'{"fallidas":>9s} {"tabla":>8s}')
estado = {}
for c in carpetas:
    gi = os.path.join(c, 'galfitm_inputs')
    ids = [d for d in os.listdir(gi) if d.isdigit() and
           os.path.exists(os.path.join(gi, d, f'{d}.feedme'))]
    salidas = [os.path.join(gi, d, f'{d}_galfitm_out.fits') for d in ids]
    hechas = [s for s in salidas if os.path.exists(s)]
    tabla = os.path.join(c, 'galfitm_resultados.csv')
    ultima = max((os.path.getmtime(s) for s in hechas), default=0)
    if not os.path.exists(tabla):
        st_tabla = 'falta'
    elif os.path.getmtime(tabla) < ultima:
        st_tabla = 'vieja'
    else:
        st_tabla = 'al dia'
    # Galaxias que ya fallaron antes (lista_fallo.txt de correr_galfitm_paralelo.sh):
    # no se reintentan cada vez, salvo con --reintentar-fallidas
    ffallo = os.path.join(c, 'lista_fallo.txt')
    fallidas = set(open(ffallo).read().split()) if os.path.exists(ffallo) else set()
    sin_salida = {d for d, sal in zip(ids, salidas) if not os.path.exists(sal)}
    pendientes = sin_salida if args.reintentar_fallidas else sin_salida - fallidas
    estado[c] = dict(n=len(ids), hechas=len(hechas), tabla=st_tabla,
                     pendientes=len(pendientes), ultima=ultima)
    print(f'{os.path.basename(c):22s} {len(ids):9d} {len(hechas):11d} {len(pendientes):11d} '
          f'{len(sin_salida & fallidas):9d} {st_tabla:>8s}')

if not args.listar and not SALTAR_PREVIOS:
    for c, e in estado.items():
        nombre = os.path.basename(c)
        if e['pendientes'] and not args.sin_ss:
            copiar_si_falta(EXE, c)
            src = os.path.join(BASE, 'correr_galfitm_paralelo.sh')
            if os.path.exists(src):                       # siempre la version de aqui
                shutil.copy2(src, os.path.join(c, 'correr_galfitm_paralelo.sh'))
            if args.reintentar_fallidas:
                os.environ['REINTENTAR'] = '1'
            correr(['bash', 'correr_galfitm_paralelo.sh', str(args.nproc)], c,
                   f'{nombre}: single-Sersic de {e["pendientes"]} galaxias')
            # la tabla solo se rehace si aparecieron salidas nuevas
            gi = os.path.join(c, 'galfitm_inputs')
            nuevas = [os.path.join(gi, d, f'{d}_galfitm_out.fits') for d in os.listdir(gi) if d.isdigit()]
            ultima = max((os.path.getmtime(x) for x in nuevas if os.path.exists(x)), default=0)
            tabla = os.path.join(c, 'galfitm_resultados.csv')
            if not os.path.exists(tabla) or os.path.getmtime(tabla) < ultima:
                e['tabla'] = 'vieja'
        if e['tabla'] != 'al dia':
            copiar_si_falta('leer_output_galfitm.py', c)
            correr([sys.executable, 'leer_output_galfitm.py'], c, f'{nombre}: leer salidas')

# =============================================================================
# 2. CORTES SOBRE TODAS LAS GALAXIAS
# =============================================================================
if not args.listar and not SALTAR_PREVIOS:
    cs = [sys.executable, 'combinar_carpetas.py'] + carpetas
    correr(cs, BASE, 'Unir tablas de todas las carpetas')
    pasos = ['procesar_muestra.py']
    if args.con_clasificacion:
        pasos += ['corregir_ext_kcor.py', 'clasificar_morfologia.py']
    for s in pasos:
        if not os.path.exists(os.path.join(BASE, s)):
            sys.exit(f'\nFalta {s} en {BASE}\nCopialo a esta carpeta (la copia dentro de '
                     f'galaxias_* no se usa).')
        correr([sys.executable, s], BASE, s)

if not os.path.exists(MUESTRA_AB):
    sys.exit('\nTodavia no hay muestra clasificada: corre sin --listar.')

# =============================================================================
# 3. QUE FALTA DE BULBO+DISCO Y STATMORPH (solo galaxias que pasan los cortes)
# =============================================================================
m = pd.read_csv(MUESTRA_AB, dtype={'objectId': str}, low_memory=False)
m = m[m.ruta_carpeta.isin(carpetas)] if 'ruta_carpeta' in m else m
ab = m[m.categoria.isin(CATEGORIAS_TRABAJO)].copy()
ab['dir'] = ab.apply(ruta_galaxia, axis=1)
ab['bd_hecho'] = [tiene_bd(d, o) for d, o in zip(ab.dir, ab.objectId)]
ab['sm_hecho'] = [tiene_sm(d, o) for d, o in zip(ab.dir, ab.objectId)]

print('\n' + '=' * 72)
print(f'Galaxias que pasan los cortes ({"+".join(CATEGORIAS_TRABAJO)}) por carpeta:')
res = ab.groupby('carpeta').agg(pasan=('objectId', 'size'),
                                bd_hecho=('bd_hecho', 'sum'),
                                sm_hecho=('sm_hecho', 'sum'))
res['bd_falta'] = res.pasan - res.bd_hecho
res['sm_falta'] = res.pasan - res.sm_hecho
print(res.to_string())
print(f'TOTAL: {len(ab)} pasan | bulbo+disco faltan {(~ab.bd_hecho).sum()} | '
      f'statmorph faltan {(~ab.sm_hecho).sum()}')

pend_bd = ab.loc[~ab.bd_hecho, 'dir'].tolist()
pend_sm = ab.loc[~ab.sm_hecho, 'dir'].tolist()
open(os.path.join(BASE, 'pendientes_bd.txt'), 'w').write('\n'.join(pend_bd) + '\n')
open(os.path.join(BASE, 'pendientes_statmorph.txt'), 'w').write('\n'.join(pend_sm) + '\n')

# Galaxias con resultados que con los cortes actuales quedaron fuera
fuera = m[~m.categoria.isin(CATEGORIAS_TRABAJO)].copy()
sobrantes = []
if len(fuera) and 'ruta_carpeta' in fuera:
    fuera['dir'] = fuera.apply(ruta_galaxia, axis=1)
    sobrantes = [(d, o) for d, o in zip(fuera.dir, fuera.objectId)
                 if os.path.isdir(d) and (tiene_bd(d, o) or tiene_sm(d, o))]
if sobrantes:
    print(f'\n{len(sobrantes)} galaxias con bulbo+disco o statmorph que ahora NO pasan los cortes')
    if args.limpiar_sobrantes and not args.listar:
        for d, o in sobrantes:
            for pat in (f'{o}_bd.feedme', f'{o}_bd.constraints', f'{o}_galfitm_bd_out.fits',
                        f'{o}_bd_galfit.*', f'galfitm_bd_{o}.log'):
                for f in glob.glob(os.path.join(d, pat)):
                    os.remove(f)
        print('  borrados sus archivos de bulbo+disco (statmorph no guarda archivos por galaxia)')
    else:
        print('  (para borrar sus archivos: --limpiar-sobrantes)')

if args.listar:
    print('\nListas escritas: pendientes_bd.txt, pendientes_statmorph.txt (no se corrio nada)')
    sys.exit(0)

# =============================================================================
# 4. CORRER LO QUE FALTA
# =============================================================================
if pend_bd and not args.sin_bd:
    if not os.path.exists(os.path.join(BASE, EXE)):
        sys.exit(f'Falta {EXE} en {BASE} (lo necesita correr_bd.sh)')
    correr([sys.executable, 'generar_bd_feedmes.py'], BASE, 'Feedmes bulbo+disco')
    # correr_bd.sh salta las que ya tienen salida: solo corre las pendientes
    correr(['bash', 'correr_bd.sh', str(args.nproc)], BASE,
           f'Bulbo+disco: {len(pend_bd)} galaxias pendientes')

if pend_sm and not args.sin_statmorph:
    if not os.path.exists(os.path.join(BASE, STATMORPH_SCRIPT)):
        sys.exit(f'Falta {STATMORPH_SCRIPT} en {BASE}')
    os.makedirs(DIR_SM, exist_ok=True)
    sello = pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')
    falta = ab[~ab.sm_hecho]
    for carpeta, grupo in falta.groupby('ruta_carpeta'):
        nombre = os.path.basename(carpeta)
        lista = os.path.join(DIR_SM, f'pendientes_{nombre}_{sello}.csv')
        grupo[['objectId']].to_csv(lista, index=False)
        out = os.path.join(DIR_SM, f'{nombre}_{sello}.csv')
        correr([sys.executable, STATMORPH_SCRIPT, carpeta, '--lista', lista,
                '--ncpu', str(args.nproc), '--out', out] + STATMORPH_ARGS,
               BASE, f'statmorph {nombre}: {len(grupo)} galaxias pendientes')

# =============================================================================
# 5. TABLA DE STATMORPH (une todas las corridas; solo galaxias que pasan los cortes)
# =============================================================================
corr = sorted(glob.glob(os.path.join(DIR_SM, '*.csv')))
corr = [f for f in corr if not os.path.basename(f).startswith('pendientes_')]
if corr:
    t = pd.concat([pd.read_csv(f, dtype={'objectId': str}, low_memory=False) for f in corr],
                  ignore_index=True)
    clave = ['objectId', 'band'] if 'band' in t else ['objectId']
    t = t.drop_duplicates(clave, keep='last')
    t = t.merge(ab[['objectId', 'carpeta']], on='objectId', how='inner')
    t.to_csv(SALIDA_SM, index=False)
    print(f'\nstatmorph: {t.objectId.nunique()} galaxias (A+B) -> {SALIDA_SM}')

print('\nListo.')

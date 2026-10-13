# -*- coding: utf-8 -*-
"""Paso 2 (clasificación visual de Claude): utilidades para trabajar por lotes y poder retomar.
  python3 paso2_claude.py orden            crea orden_paso2_claude.txt (una vez)
  python3 paso2_claude.py lote [N]         muestra las próximas N sin clasificar, con las cantidades de las reglas
  python3 paso2_claude.py guardar <txt>    agrega un lote: líneas  objectId|visual_ok|letras de etiquetas|nota
  python3 paso2_claude.py estado
Nunca lee las columnas visual* de inspeccion_visual.csv (solo objectId) ni escribe en ese archivo."""
import sys, os
from datetime import datetime
import numpy as np, pandas as pd
D = os.path.dirname(os.path.abspath(__file__))
AUX, ORD = os.path.join(D, 'reglas_aux_claude.csv'), os.path.join(D, 'orden_paso2_claude.txt')
OUT, USR, PNG = os.path.join(D, 'inspeccion_visual_claude.csv'), os.path.join(D, 'inspeccion_visual.csv'), os.path.join(D, 'inspeccion_png')
NOMBRE = {1: 'aceptada', 0: 'rechazada', 2: 'dudosa', 3: 'una_componente'}
ETQ = {'b': 'barra', 'o': 'anillo', 'v': 'vecino', 's': 'estructura', 'e': 'artefacto', 'i': 'intercambiadas'}
COLS = ['objectId', 'visual', 'visual_ok', 'visual_etiquetas', 'categoria_bd', 'banda_vista', 'fecha', 'clasificador', 'nota']
aux = pd.read_csv(AUX, dtype={'objectId': str}).set_index('objectId')
usuaria = set(pd.read_csv(USR, dtype={'objectId': str}, usecols=['objectId']).objectId)


# Trabajo en paralelo: con PARTE=k/K cada proceso toma orden[k::K] y guarda en paso2_partes/parte_k.csv
# (así no se pisan); "unir" junta las partes en inspeccion_visual_claude.csv.
PARTES = os.path.join(D, 'paso2_partes')
PARTE = os.environ.get('PARTE')
K_, N_ = (int(x) for x in PARTE.split('/')) if PARTE else (None, None)
DEST = os.path.join(PARTES, f'parte_{K_}.csv') if PARTE else OUT


def leer_csv(p):
    if os.path.exists(p):
        return pd.read_csv(p, dtype={'objectId': str, 'visual_etiquetas': str, 'nota': str})
    return pd.DataFrame(columns=COLS)


def hechas():
    import glob
    return pd.concat([leer_csv(p) for p in [OUT] + sorted(glob.glob(os.path.join(PARTES, 'parte_*.csv')))],
                     ignore_index=True).drop_duplicates('objectId', keep='last')


def pista(r):
    if r.categoria_bd == 'disco_puro':
        return 'DP(3/0)'
    p = []
    x = r.r_cruce / r.fwhm_arcsec
    if 0.9 <= x <= 1.1:
        p.append('2?fig')
    elif x < 0.9:
        p.append('2S' if r.BT < 0.2 else ('2A' if r.ratio_re < 0.5 else '2V'))
    if r.n_b < 1.5 and r.ratio_re < 0.5:
        p.append('3a')
    if r.n_b < 0.65 and r.ratio_re >= 0.5 and r.BT > 0.5:
        p.append('3s')
    if r.BT < 0.11 and r.q_b < 0.4:
        p.append('4n')
    return ' '.join(p)


cmd = sys.argv[1]
if cmd == 'orden':
    if os.path.exists(ORD):
        sys.exit('ya existe ' + ORD)
    s = aux[aux.categoria_bd.isin(['bd_ok', 'bulbo_no_resuelto', 'disco_puro']) & (aux.una_comp_motivo.fillna('') == '')]
    s = s[~s.index.isin(usuaria)].sample(frac=1, random_state=42)
    open(ORD, 'w').write('\n'.join(s.index) + '\n')
    print(len(s), s.categoria_bd.value_counts().to_dict())
    sys.exit()
orden = [l.strip() for l in open(ORD) if l.strip()]
mios = orden[K_::N_] if PARTE else orden
h = hechas()
if cmd == 'estado':
    print(f'{len(h)} de {len(orden)}');  print(h.visual.value_counts().to_string() if len(h) else '')
elif cmd == 'lote':
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    ya = set(h.objectId)
    pend = [o for o in mios if o not in ya and o not in usuaria]
    n_pend = len(pend)
    pend = ([o for o in pend if os.path.exists(os.path.join(PNG, o + '_insp.png'))] or pend)[:n]   # primero las que ya tienen figura
    print(f'hechas en total {len(ya)} / {len(orden)}; pendientes de esta parte: {n_pend}')
    for o in pend:
        r = aux.loc[o]
        png = '' if os.path.exists(os.path.join(PNG, o + '_insp.png')) else '  SIN PNG'
        print(f'{o} {r.categoria_bd[:5]} BT={r.BT:.2f} n={r.n_b:.1f} re={r.ratio_re:.2f} qb={r.q_b:.2f} qd={r.q_d:.2f} dPA={r.dPA:.0f} x={r.r_cruce:.2f}/{r.fwhm_arcsec:.2f} [{pista(r)}]{png}')
elif cmd == 'guardar':
    filas = []
    for l in open(sys.argv[2]):
        if not l.strip():
            continue
        o, ok, et, nota = [x.strip() for x in l.rstrip('\n').split('|', 3)]
        ok = int(ok)
        assert o in orden and o not in usuaria, f'{o}: no es del paso 2'
        assert ok in NOMBRE and all(e in ETQ for e in et), l
        assert not (aux.loc[o].categoria_bd == 'disco_puro' and ok == 1), f'{o}: disco_puro solo admite 3 o 0'
        filas.append([o, NOMBRE[ok], ok, ';'.join(sorted(ETQ[e] for e in set(et))), aux.loc[o].categoria_bd, 'r',
                      datetime.now().isoformat(timespec='seconds'), 'claude', nota])
    n = pd.DataFrame(filas, columns=COLS)
    assert n.objectId.is_unique
    os.makedirs(PARTES, exist_ok=True)
    d = leer_csv(DEST)
    t = pd.concat([d[~d.objectId.isin(n.objectId)], n], ignore_index=True)
    t.to_csv(DEST + '.tmp', index=False); os.replace(DEST + '.tmp', DEST)
    print(f'guardadas {len(n)} en {os.path.basename(DEST)} (lleva {len(t)})')
elif cmd == 'unir':
    if os.path.exists(OUT) and len(leer_csv(OUT)) >= len(orden) and '--forzar' not in sys.argv:
        sys.exit('inspeccion_visual_claude.csv ya está completo y tiene correcciones posteriores (reclasificaciones): '
                 'no se vuelve a unir desde paso2_partes/ (usa --forzar si de verdad quieres pisarlo)')
    t = h.set_index('objectId').loc[[o for o in orden if o in set(h.objectId)]].reset_index()[COLS]
    t.to_csv(OUT + '.tmp', index=False); os.replace(OUT + '.tmp', OUT)
    print(f'{OUT}: {len(t)} de {len(orden)}');  print(t.visual.value_counts().to_string())

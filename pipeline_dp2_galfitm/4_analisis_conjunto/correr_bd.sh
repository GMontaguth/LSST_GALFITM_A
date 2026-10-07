#!/bin/bash
# Corre GalfitM bulbo + disco en las carpetas de lista_bd.txt
# (generada por generar_bd_feedmes.py). Cada ajuste corre DENTRO de la
# carpeta de su galaxia, asi los archivos galfit.NN no se mezclan y se
# puede correr en paralelo.
#
# Uso:  ./correr_bd.sh          # 1 galaxia a la vez
#       ./correr_bd.sh 4        # 4 en paralelo
# Es reanudable: las galaxias que ya tienen salida se saltan.

BASE="$(cd "$(dirname "$0")" && pwd)"
export GALFITM="$BASE/galfitm-1.4.4-linux-x86_64"
export SUFIJO="_bd"
NPROC="${1:-1}"
LISTA="$BASE/lista_bd.txt"

[ -x "$GALFITM" ] || { echo "No encuentro o no es ejecutable: $GALFITM"; exit 1; }
[ -f "$LISTA" ]   || { echo "No existe $LISTA (corre generar_bd_feedmes.py)"; exit 1; }

correr_una() {
    dir="$1"
    id=$(basename "$dir")
    cd "$dir" || { echo "FALLO $id (no existe la carpeta)"; return; }
    salida="${id}_galfitm${SUFIJO}_out.fits"
    if [ -f "$salida" ]; then
        echo "SALTADA $id (ya tiene salida)"
        return
    fi
    marca=$(mktemp -p . .inicio_XXXX)
    "$GALFITM" "${id}${SUFIJO}.feedme" > "galfitm${SUFIJO}_${id}.log" 2>&1
    # parametros finales: galfit.NN creados en ESTA corrida -> <ID>_bd_galfit.NN
    for f in $(find . -maxdepth 1 -name 'galfit.[0-9]*' -newer "$marca" -printf '%f\n'); do
        mv "$f" "${id}${SUFIJO}_${f}"
    done
    rm -f "$marca"
    if [ -f "$salida" ]; then echo "OK $id"; else echo "FALLO $id (ver galfitm${SUFIJO}_${id}.log)"; fi
}
export -f correr_una

total=$(grep -c . "$LISTA")
echo "Corriendo $total galaxias con $NPROC proceso(s)..."
grep . "$LISTA" | xargs -P "$NPROC" -I{} bash -c 'correr_una "$@"' _ {} | tee "$BASE/progreso_bd.txt"

echo ""
echo "===== TERMINADO ====="
echo "OK:       $(grep -c '^OK'      "$BASE/progreso_bd.txt")"
echo "Saltadas: $(grep -c '^SALTADA' "$BASE/progreso_bd.txt")"
echo "Fallidas: $(grep -c '^FALLO'   "$BASE/progreso_bd.txt")"
grep '^FALLO' "$BASE/progreso_bd.txt" > "$BASE/fallidas_bd.txt"

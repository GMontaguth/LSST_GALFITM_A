#!/bin/bash
# Corre GalfitM (single-Sersic) en todas las galaxias de galfitm_inputs/
# con varios procesos a la vez.
#
# Uso (desde galaxias_X_Y, junto al ejecutable):
#   ./correr_galfitm_paralelo.sh        # 1 proceso
#   ./correr_galfitm_paralelo.sh 3      # 3 procesos en paralelo
#
# Por que no basta con lanzar 3 veces correr_galfitm.sh: los feedme usan rutas
# 'galfitm_inputs/<ID>/...' y GalfitM escribe galfit.01, galfit.02... en la
# carpeta desde donde se ejecuta. Con varios procesos en la misma carpeta esos
# archivos se pisan. Aqui cada ajuste corre en su propia carpeta temporal, con
# un enlace a galfitm_inputs/, asi las rutas siguen valiendo y nada se mezcla.
#
# Es reanudable: las galaxias que ya tienen <ID>_galfitm_out.fits se saltan.

BASE="$(cd "$(dirname "$0")" && pwd)"
export GALFITM="$BASE/galfitm-1.4.4-linux-x86_64"
export INPUTS="$BASE/galfitm_inputs"
NPROC="${1:-1}"

[ -x "$GALFITM" ] || { echo "No encuentro o no es ejecutable: $GALFITM"; exit 1; }
[ -d "$INPUTS" ]  || { echo "No existe $INPUTS"; exit 1; }

correr_una() {
    id="$1"
    gdir="$INPUTS/$id"
    feedme="galfitm_inputs/$id/$id.feedme"
    salida="$gdir/${id}_galfitm_out.fits"

    [ -f "$gdir/$id.feedme" ] || { echo "SIN_FEEDME $id"; return; }
    [ -f "$salida" ] && { echo "SALTADA $id"; return; }

    tmp=$(mktemp -d)                         # carpeta propia de este ajuste
    ln -s "$INPUTS" "$tmp/galfitm_inputs"
    ( cd "$tmp" && "$GALFITM" "$feedme" > "$gdir/galfitm_$id.log" 2>&1 )

    # parametros finales -> carpeta de la galaxia como <ID>_galfit.NN
    for f in "$tmp"/galfit.[0-9]*; do
        [ -e "$f" ] && mv "$f" "$gdir/${id}_$(basename "$f")"
    done
    [ -f "$tmp/fit.log" ] && cat "$tmp/fit.log" >> "$gdir/fit.log"
    rm -rf "$tmp"

    if [ -f "$salida" ]; then echo "OK $id"; else echo "FALLO $id"; fi
}
export -f correr_una

ls "$INPUTS" | grep -E '^[0-9]+$' > "$BASE/.lista_ss.txt"
# Las que ya fallaron antes no se reintentan (REINTENTAR=1 para reintentarlas)
FALLO_PREV="$BASE/lista_fallo.txt"
touch "$FALLO_PREV"
cut -d' ' -f1 "$FALLO_PREV" | grep -E '^[0-9]+$' | sort -u > "$BASE/.fallo_prev.txt"
if [ "${REINTENTAR:-0}" != "1" ] && [ -s "$BASE/.fallo_prev.txt" ]; then
    grep -vxF -f "$BASE/.fallo_prev.txt" "$BASE/.lista_ss.txt" > "$BASE/.lista_tmp.txt"
    echo "Se omiten $(wc -l < "$BASE/.fallo_prev.txt") galaxias que ya fallaron (REINTENTAR=1 para reintentarlas)"
    mv "$BASE/.lista_tmp.txt" "$BASE/.lista_ss.txt"
fi
total=$(wc -l < "$BASE/.lista_ss.txt")
echo "Corriendo $total galaxias con $NPROC proceso(s)..."

xargs -P "$NPROC" -I{} bash -c 'correr_una "$@"' _ {} < "$BASE/.lista_ss.txt" \
    | tee "$BASE/progreso_ss.txt"
rm -f "$BASE/.lista_ss.txt"

echo ""
echo "===== TERMINADO ====="
echo "OK:        $(grep -c '^OK'         "$BASE/progreso_ss.txt")"
echo "Saltadas:  $(grep -c '^SALTADA'    "$BASE/progreso_ss.txt")"
echo "Fallidas:  $(grep -c '^FALLO'      "$BASE/progreso_ss.txt")"
echo "Sin feedme:$(grep -c '^SIN_FEEDME' "$BASE/progreso_ss.txt")"
# lista_fallo.txt = fallidas de esta corrida (+ las anteriores que no se reintentaron)
grep '^FALLO' "$BASE/progreso_ss.txt" | cut -d' ' -f2 > "$BASE/.fallo_nuevo.txt"
if [ "${REINTENTAR:-0}" = "1" ]; then
    sort -u "$BASE/.fallo_nuevo.txt" > "$FALLO_PREV"
else
    sort -u "$BASE/.fallo_prev.txt" "$BASE/.fallo_nuevo.txt" > "$FALLO_PREV"
fi
rm -f "$BASE/.fallo_prev.txt" "$BASE/.fallo_nuevo.txt"

#!/bin/bash
# Junta las carpetas de galaxias de un rango de lotes en una carpeta de trabajo,
# con la misma estructura que galaxias_1_9:
#
#   galaxias_10_22/
#       galfitm_inputs/<objectId>/...     <- carpetas de las galaxias
#       galfitm-1.4.4-linux-x86_64        <- copiado
#       *.py, *.sh                        <- scripts copiados
#
# Uso (desde la carpeta que contiene los galfitm_dp2_z01_loteXX):
#   ./crear_carpeta_lotes.sh 10 22            # solo muestra lo que haria
#   ./crear_carpeta_lotes.sh 10 22 --mover    # lo hace (MUEVE las carpetas)
#
# Busca las galaxias hasta 3 niveles dentro de cada lote (algunos tar se
# descomprimen con una subcarpeta extra) y solo toma carpetas con nombre
# numerico que tengan un .feedme. Si un objectId aparece en dos lotes, se
# conserva el primero y se avisa.

INI="$1"; FIN="$2"; MODO="$3"
[ -z "$INI" ] || [ -z "$FIN" ] && { echo "Uso: $0 <lote_ini> <lote_fin> [--mover]"; exit 1; }

DEST="galaxias_${INI}_${FIN}"
ORIGEN_SCRIPTS="galaxias_1_9"          # de aqui se copian el ejecutable y los scripts
EXE="galfitm-1.4.4-linux-x86_64"

lista=$(mktemp)
for n in $(seq "$INI" "$FIN"); do
    nn=$(printf '%02d' "$n")
    for lote in galfitm_dp2_z01_lote${nn} galfitm_dp2_z01_lote${nn}.*; do
        [ -d "$lote" ] || continue
        find "$lote" -mindepth 1 -maxdepth 3 -type d -regextype posix-extended \
             -regex '.*/[0-9]+$' -exec sh -c 'ls "$1"/*.feedme >/dev/null 2>&1 && echo "$1"' _ {} \;
    done
done > "$lista"

total=$(wc -l < "$lista")
unicos=$(xargs -n1 basename < "$lista" | sort -u | wc -l)
echo "Lotes $INI a $FIN"
echo "  carpetas de galaxias encontradas : $total"
echo "  objectId distintos               : $unicos"
[ "$total" -ne "$unicos" ] && echo "  AVISO: $((total - unicos)) repetidas entre lotes (se conserva la primera)"
echo "  por lote:"
sed -E 's#^(galfitm_dp2_z01_lote[0-9.]+)/.*#\1#' "$lista" | sort | uniq -c | sed 's/^/    /'

if [ "$MODO" != "--mover" ]; then
    echo ""
    echo "No se movio nada. Para hacerlo:  $0 $INI $FIN --mover"
    rm -f "$lista"; exit 0
fi

mkdir -p "$DEST/galfitm_inputs"
movidas=0
while read -r d; do
    id=$(basename "$d")
    if [ -e "$DEST/galfitm_inputs/$id" ]; then
        echo "  repetida, se deja en su lote: $d"
    else
        mv "$d" "$DEST/galfitm_inputs/" && movidas=$((movidas + 1))
    fi
done < "$lista"
rm -f "$lista"

# ejecutable y scripts
if [ -d "$ORIGEN_SCRIPTS" ]; then
    [ -e "$DEST/$EXE" ] || cp "$ORIGEN_SCRIPTS/$EXE" "$DEST/" 2>/dev/null
    chmod +x "$DEST/$EXE" 2>/dev/null
    for f in correr_galfitm.sh leer_output_galfitm.py filtrar_gri_local.py procesar_muestra.py \
             corregir_ext_kcor.py calc_kcor_fun.py clasificar_morfologia.py \
             seleccionar_discos.py generar_bd_feedmes.py correr_bd.sh; do
        [ -f "$ORIGEN_SCRIPTS/$f" ] && [ ! -e "$DEST/$f" ] && cp "$ORIGEN_SCRIPTS/$f" "$DEST/"
    done
    chmod +x "$DEST"/*.sh 2>/dev/null
fi

echo ""
echo "Listo: $movidas galaxias en $DEST/galfitm_inputs/"
echo "Archivos copiados a $DEST/: $(ls "$DEST" | grep -v galfitm_inputs | tr '\n' ' ')"

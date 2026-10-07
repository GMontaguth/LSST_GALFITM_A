# Pipeline GalfitM + statmorph + bulbo/disco — LSST DP2 (z < 0.1)

Desde la descarga de las imágenes en el RSP (paso 0) hasta la medición no
paramétrica con **statmorph** y la descomposición **bulbo + disco** con GalfitM,
para galaxias con redshift espectroscópico a z < 0.1.

El pipeline está pensado para trabajar **por partes**: las galaxias se descargan
por lotes, se agrupan en carpetas de trabajo (`galaxias_1_9`, `galaxias_10_22`...)
y cada vez que se suma una carpeta, un solo comando (`actualizar_todo.py`) corre
**solo lo que falta**.

---

## Resumen

| Paso | Dónde | Script | Qué hace |
|---|---|---|---|
| 0 | RSP (notebook) | `z01_dp2_galfitm_lotes_v8.ipynb` | Descarga stamps, PSF, sigma, máscara y feedme por lotes (solo galaxias con g, r e i) |
| 1 | `Galfitm/` | `crear_carpeta_lotes.sh` | Junta varios lotes en una carpeta de trabajo `galaxias_X_Y` |
| 2 | `galaxias_X_Y/` | `filtrar_gri_local.py` | (solo lotes bajados con la v7) quita galaxias sin g, r, i útiles |
| 3 | `galaxias_X_Y/` | `correr_galfitm_paralelo.sh` | Ajuste single-Sérsic multibanda con GalfitM |
| 4 | `galaxias_X_Y/` | `leer_output_galfitm.py` | Lee las salidas → `galfitm_resultados.csv` |
| 5 | `analisis_conjunto/` | `combinar_carpetas.py` | Une las tablas de todas las carpetas (columna `carpeta`) |
| 6 | `analisis_conjunto/` | `procesar_muestra.py` | Cortes de calidad, cruce con LSST y categorías A/B/C/D |
| 7 | `analisis_conjunto/` | `statmorph_completo.py` | statmorph-lsst para las galaxias A + B |
| 8 | `analisis_conjunto/` | `generar_bd_feedmes.py` + `correr_bd.sh` | Bulbo + disco con GalfitM para las galaxias A + B |

Los pasos 3 a 8 los coordina **`actualizar_todo.py`** (ver [Uso diario](#uso-diario-actualizar_todopy)).

---

## Estructura de carpetas

```
Galfitm/
├── galfitm_dp2_z01_loteXX.tar.gz        # lotes descargados del RSP
├── crear_carpeta_lotes.sh
├── statmorph-lsst/                      # librería statmorph-lsst (clonada)
├── galaxias_1_9/                        # carpeta de trabajo (lotes 1-9)
│   ├── galfitm_inputs/<objectId>/       # una carpeta por galaxia
│   ├── galfitm-1.4.4-linux-x86_64
│   ├── correr_galfitm_paralelo.sh       # (se copia sola)
│   ├── leer_output_galfitm.py           # (se copia sola)
│   └── galfitm_resultados.csv
├── galaxias_10_22/                      # carpeta de trabajo (lotes 10-22)
│   └── ...
└── analisis_conjunto/                   # análisis de TODAS las carpetas
    ├── actualizar_todo.py
    ├── combinar_carpetas.py
    ├── procesar_muestra.py
    ├── statmorph_completo.py
    ├── generar_bd_feedmes.py
    ├── correr_bd.sh
    ├── correr_galfitm_paralelo.sh
    ├── leer_output_galfitm.py
    ├── galfitm-1.4.4-linux-x86_64
    └── analisis_galfitm/                # resultados (tablas y figuras)
```

Contenido de cada carpeta de galaxia (`galfitm_inputs/<objectId>/`):

| Archivo | Origen |
|---|---|
| `<ID>_sci_<b>.fits`, `<ID>_sig_<b>.fits`, `<ID>_psf_<b>.fits`, `<ID>_mask.fits` | paso 0 |
| `<ID>.feedme` | paso 0 (single-Sérsic) |
| `<ID>_galfitm_out.fits`, `<ID>_galfit.01`, `galfitm_<ID>.log` | paso 3 |
| `<ID>_bd.feedme`, `<ID>_bd.constraints` | paso 8 |
| `<ID>_galfitm_bd_out.fits`, `<ID>_bd_galfit.01`, `galfitm_bd_<ID>.log` | paso 8 |

## Requisitos

- Python ≥ 3.9 con `numpy`, `pandas`, `matplotlib`, `scipy`, `astropy`, `pyarrow`
- `photutils` y la librería **statmorph-lsst** (`pip install -e Galfitm/statmorph-lsst`,
  o dejarla en `Galfitm/statmorph-lsst`: `statmorph_completo.py` la encuentra sola)
- GalfitM 1.4.4 (`galfitm-1.4.4-linux-x86_64`)

---

## Paso 0 — Descarga en el RSP (`z01_dp2_galfitm_lotes_v8.ipynb`)

Notebook para el Rubin Science Platform. Para cada galaxia del catálogo
espectroscópico cruzado con LSST:

- recorta un stamp por banda (imagen, varianza → sigma), la PSF y genera la
  máscara de vecinos (SEP);
- escribe el feedme single-Sérsic con valores iniciales de LSST;
- trabaja por **lotes de 800 galaxias** (`LOTE = 1, 2, 3...`) y empaqueta cada
  lote en `galfitm_dp2_z01_loteXX.tar.gz`.

**Filtro g, r, i (v8).** Solo se descargan galaxias con las tres bandas
(`BANDAS_OBLIGATORIAS = ['g', 'r', 'i']`): se piden primero sus PSF y, si falta
alguna, la galaxia se descarta sin bajar el resto. Si el stamp de g, r o i tiene
más del 30 % de píxeles vacíos (`FRAC_VACIA_MAX`), la galaxia también se descarta.
Las bandas u, z, y se descargan cuando existen, pero no son obligatorias.

El filtro está dentro del bucle (no en el catálogo), así que **la numeración de
los lotes no cambia** entre versiones: los lotes 1–15 se bajaron con la v7 y del
16 en adelante con la v8.

> En la última celda, `CONFIRMAR_BORRADO` borra el lote del servidor. Conviene
> dejarlo en `False` y activarlo a mano solo después de descargar el tar.

## Paso 1 — Agrupar lotes (`crear_carpeta_lotes.sh`)

```bash
cd Galfitm
./crear_carpeta_lotes.sh 10 22            # solo muestra lo que haría
./crear_carpeta_lotes.sh 10 22 --mover    # crea galaxias_10_22/ y mueve las galaxias
```

Busca las carpetas de galaxias (nombre numérico y con `.feedme`) dentro de los
lotes indicados, hasta 3 niveles de profundidad (algunos tar traen una subcarpeta
extra), las **mueve** a `galaxias_X_Y/galfitm_inputs/` y copia el ejecutable y los
scripts desde `galaxias_1_9`. Si un objectId aparece en dos lotes, conserva el
primero y avisa.

## Paso 2 — Filtro g, r, i local (`filtrar_gri_local.py`)

Solo para lotes descargados con la v7 (1–15). Aplica el mismo filtro que la v8:

```bash
cd galaxias_10_22
python3 filtrar_gri_local.py            # crea sin_gri.txt (no borra)
python3 filtrar_gri_local.py --borrar   # borra esas carpetas
```

No es obligatorio (el paso 6 descarta igual esas galaxias), pero libera espacio
y evita ajustarlas.

## Paso 3 — Single-Sérsic (`correr_galfitm_paralelo.sh`)

```bash
cd galaxias_10_22
./correr_galfitm_paralelo.sh 3       # 3 procesos en paralelo
```

- Cada ajuste corre en una **carpeta temporal propia** con un enlace a
  `galfitm_inputs/`: las rutas del feedme (`galfitm_inputs/<ID>/...`) siguen
  valiendo y los `galfit.NN` de distintos procesos no se mezclan.
- Los parámetros finales se guardan como `<ID>_galfit.01` en la carpeta de la galaxia.
- **Reanudable**: salta las galaxias con `<ID>_galfitm_out.fits`.
- Las que fallan van a `lista_fallo.txt` y **no se reintentan** en corridas
  siguientes (`REINTENTAR=1 ./correr_galfitm_paralelo.sh 3` para reintentarlas).

Receta del feedme: un Sérsic más cielo, con polinomios de Chebyshev en λ
(posición constante; R_e, n, b/a y PA lineales; magnitud libre por banda),
zeropoint 31.4 (flujos en nJy → AB) y escala 0.2″/px.

## Paso 4 — Leer salidas (`leer_output_galfitm.py`)

```bash
python3 leer_output_galfitm.py
```

Lee el header de cada `<ID>_galfitm_out.fits` y crea `galfitm_resultados.csv`
(una fila por galaxia): `MAG_b, RE_b, N_b, AR_b, PA_b, XC_b, YC_b` y sus errores
por banda, `SKY_b`, `chi2nu`, `flag_b` (parámetros marcados con `*` por GalfitM)
y `estado` (`ok`, `sin_output`, `error`). Hace además figuras input/modelo/residuo
de 100 galaxias al azar (`galfitm_imagenes/`).

---

## Paso 5 — Unir carpetas (`combinar_carpetas.py`)

Une los `galfitm_resultados.csv` de todas las carpetas `galaxias_*` y agrega:

- `carpeta`: nombre de la carpeta de origen;
- `ruta_carpeta`: ruta completa (los scripts la usan para encontrar los FITS).

Si un objectId está en dos carpetas, se conserva el ajuste `ok` con menor χ²/ν.

**Los cortes se hacen sobre la muestra completa**, no por carpeta, porque varios
umbrales se calculan con la población (dispersión de colores, offsets con LSST):
por carpeta saldrían umbrales distintos para galaxias equivalentes.

## Paso 6 — Cortes y muestra confiable (`procesar_muestra.py`)

Tres etapas en un script (`python3 procesar_muestra.py [2] [3] [4]`).

### 6a. Cortes de calidad de GalfitM

Por galaxia: `0.5 ≤ χ²/ν ≤ 5`.

Por banda (solo esa banda pasa a NaN):

| Motivo | Condición |
|---|---|
| `flag_galfitm` | GalfitM marcó un parámetro con `*` |
| `sin_error` | error de magnitud vacío o ≤ 0 |
| `mag_fuera` | magnitud fuera de [10, 25] |
| `re_fuera` | R_e fuera de [0.5, 150] px |
| `n_limite` | n fuera de [0.21, 7.9] |
| `ar_bajo` | b/a < 0.1 |
| `color_outlier` | magnitud > 5σ (y > 1 mag) fuera del color típico de la población |
| `imagen_vacia` | > 30 % de píxeles vacíos en la imagen INPUT |
| `fondo_malo` | variación de fondo a gran escala > 3 × ruido |

Las métricas de imagen se guardan en `calidad_imagenes.csv` y **solo se miden
las galaxias nuevas** en cada corrida.

### 6b. Cruce con LSST

Por `objectId` con el catálogo (`LSST_PARQUET`). Si hay varios espectros por
objeto, se conserva `f_z = 1` y el de menor distancia. Se convierte
`{b}_sersicFlux` (nJy) a magnitud AB con el mismo zeropoint (31.4).

### 6c. Categorías A / B / C / D

**LSST no se toma como la verdad.** MultiProFit es otro ajuste paramétrico
(un solo n para todas las bandas, n en [0.5, 6], afectado por el *deblender*).
Se usa como segunda medida y solo para detectar diferencias **enormes**:

| Parámetro | Desacuerdo si la diferencia (respecto al offset mediano) supera |
|---|---|
| Magnitud | 1.0 mag (y ≥ 5σ) |
| R_e | 0.30 dex = factor 2 (y ≥ 5σ) |
| n | 0.30 dex = factor 2 (y ≥ 5σ); no se compara si n de LSST está en su límite |
| b/a | 0.25 (y ≥ 5σ) |

Cuando hay desacuerdo, decide **la imagen** (independiente de ambos modelos),
en una elipse de 2 R_e: flujo de apertura sobre los píxeles, cociente flujo
modelo / flujo imagen y *residual flux fraction* (RFF, Hoyos et al. 2011).

| Categoría | Condición |
|---|---|
| **A** | Coincide con LSST |
| **B** | Desacuerdo enorme, pero GalfitM reproduce la imagen (RFF < 0.10 y flujo modelo/imagen en 0.9–1.1) |
| **C** | Desacuerdo enorme y GalfitM tampoco reproduce la imagen |
| **D** | No pasa los cortes de 6a, le falta g, r o i buena, tiene < 3 bandas buenas, o espectro no confiable (`f_z ≠ 1`, clase ≠ GALAXY, cruce > 1″) |

**Las galaxias A + B son las que pasan a statmorph y bulbo + disco.**

Además se mide la curva de crecimiento en la imagen para comprobar qué R_e
encierra la mitad de la luz (figura `21_radio_efectivo.png`).

Salidas en `analisis_galfitm/muestra/`: `muestra_clasificada.csv` (todas, con
`categoria`, `motivo_categoria`, `problemas_lsst`), `muestra_AB.csv`,
`metricas_imagen.csv` (caché) y figuras 12–14 y 21.

---

## Paso 7 — statmorph (`statmorph_completo.py`)

Medición no paramétrica con **statmorph-lsst** (CAS, Gini–M20, MID, radios de
Petrosian, Sérsic) en g, r, i, usando la sigma como mapa de pesos, la PSF y la
máscara de cada galaxia.

```bash
python3 statmorph_completo.py <carpeta_galaxias> --lista pendientes.csv \
        --ncpu 3 --out salida.csv --sin-doble --png 0
```

| Opción | Uso |
|---|---|
| `--lista` | CSV con los `objectId` a medir (columna `--col`, por defecto `objectId`) |
| `--sin-doble` | sin el ajuste Sérsic doble (más rápido) |
| `--png N` | figuras de control de las N primeras galaxias (0 = ninguna) |

Guarda un parcial cada 10 galaxias. Normalmente **no se corre a mano**:
`actualizar_todo.py` arma la lista de pendientes de cada carpeta y lo lanza.

Resultados:
- cada corrida → `analisis_galfitm/statmorph_corridas/<carpeta>_<fecha>.csv`;
- tabla unida (solo A + B, sin duplicados, con `carpeta`) →
  **`analisis_galfitm/statmorph_resultados.csv`**.

Una galaxia cuenta como hecha si su `objectId` aparece en cualquier CSV de
`statmorph_corridas/` (para aprovechar corridas anteriores, copiarlas ahí).

## Paso 8 — Bulbo + disco

### Feedmes (`generar_bd_feedmes.py`)

Para cada galaxia A + B, a partir de su feedme single-Sérsic. Receta de
**Vika et al. (2014, MNRAS 444, 3603, sec. 2.2.2)** (`INICIO = 'vika2014'`):

| | Bulbo (Sérsic) | Disco (exponencial) |
|---|---|---|
| Magnitud inicial | m_ss + 0.75 | m_ss + 0.65 |
| Radio inicial | R_e = 0.5 R_e,ss | R_e = R_e,ss (R_s = R_e / 1.678) |
| n inicial | n_ss | 1 (fijo) |
| b/a, PA iniciales | 0.8, 10° | los del single-Sérsic |
| Variación con λ | R_e, n, b/a, PA constantes; magnitud libre por banda | ídem |

- Mismo centro para ambas componentes; cielo fijo.
- Restricciones (`<ID>_bd.constraints`): magnitudes 5–35, R_e 0.04–600″,
  n 0.1–15, centro común con desplazamiento ≤ √(s/8) px.
- Las bandas que no pasaron los cortes **se quitan del ajuste**.
- Exige g, r, i y al menos 3 bandas buenas.
- `NO_SOBRESCRIBIR = True`: no toca feedmes existentes (protege ajustes en curso).
- Rutas relativas a la carpeta de la galaxia (cada ajuste corre ahí dentro).

Recetas alternativas: `'v2'` (versión DP1) y `'proporcional'`.

### Ajuste (`correr_bd.sh`)

```bash
./correr_bd.sh 3        # 3 galaxias en paralelo
```

Lee `lista_bd.txt`, corre GalfitM **dentro de la carpeta de cada galaxia**,
renombra `galfit.NN` → `<ID>_bd_galfit.NN`, salta las que ya tienen
`<ID>_galfitm_bd_out.fits` (reanudable) y deja `progreso_bd.txt` y
`fallidas_bd.txt`.

---

## Uso diario: `actualizar_todo.py`

Desde `analisis_conjunto/`. Recorre las carpetas `galaxias_*` (junto a esta
carpeta o un nivel arriba) y corre solo lo pendiente:

1. single-Sérsic de galaxias sin salida (no reintenta las fallidas) y
   `leer_output_galfitm.py` si hay salidas nuevas;
2. `combinar_carpetas.py` y `procesar_muestra.py` (cachés: solo mide lo nuevo);
3. statmorph y bulbo + disco **solo** para las A + B que no los tienen;
4. une los resultados de statmorph.

```bash
python3 actualizar_todo.py --listar                 # qué falta, sin correr nada
python3 actualizar_todo.py --nproc 3                # todo lo pendiente
python3 actualizar_todo.py --solo-statmorph         # solo statmorph (usa los cortes ya hechos)
python3 actualizar_todo.py --sin-bd                 # todo menos bulbo + disco
python3 actualizar_todo.py --sin-statmorph          # todo menos statmorph
python3 actualizar_todo.py --sin-ss                 # no lanzar single-Sérsic
python3 actualizar_todo.py --reintentar-fallidas    # reintentar single-Sérsic fallidos
python3 actualizar_todo.py --limpiar-sobrantes      # borrar bulbo+disco de galaxias que salieron de A+B
```

`--listar` muestra, por carpeta, cuántas galaxias pasan los cortes y cuántas
tienen o les faltan bulbo + disco (`bd_hecho`, `bd_falta`) y statmorph
(`sm_hecho`, `sm_falta`).

Para sumar una carpeta nueva basta con crearla (paso 1) y volver a correr
`actualizar_todo.py`.

---

## Buenas prácticas al correr

- **Ver el avance en pantalla** (sin `nohup`):
  ```bash
  OMP_NUM_THREADS=1 PYTHONUNBUFFERED=1 python3 actualizar_todo.py --solo-statmorph --nproc 3 2>&1 | tee statmorph_log.txt
  ```
  No cerrar esa terminal mientras corre.
- **En segundo plano**: `nohup ... > log.txt 2>&1 &` y `tail -f log.txt`.
- **Detener: Ctrl+C. Nunca Ctrl+Z** (pausa el proceso; queda en estado `T`
  ocupando memoria). Todos los pasos son reanudables.
- **No lanzar dos veces lo mismo.** Antes de lanzar, revisar:
  ```bash
  ps -eo pid,stat,lstart,args | grep -E "statmorph_completo|correr_bd|galfitm-1.4.4" | grep -v grep
  ```
  `R`/`S` = activo, `T` = pausado (eliminar con `kill -9 PID`), `Z` = zombi (inofensivo).
- **No correr bulbo + disco desde dos lugares a la vez** (usar `--sin-bd` si ya
  hay un `correr_bd.sh` activo).
- `OMP_NUM_THREADS=1` evita que statmorph compita por núcleos con GalfitM.
  Con `nproc` núcleos, no superar `nproc − 1` procesos en total.

---

## Resultados de validación (muestra actual)

Con 7466 galaxias (`galaxias_1_9` + `galaxias_10_22`):

| Categoría | N |
|---|---|
| A | 3939 |
| B | 900 |
| C | 383 |
| D | 2244 (1334 por faltar g, r o i buena) |

- En los desacuerdos, la magnitud de LSST es más débil que el flujo de apertura
  medido en la imagen en el **75 %** de los casos (GalfitM: 1.5 %): cuando difieren,
  quien pierde flujo suele ser LSST (galaxias con `shape_flag` o fragmentadas).
- **Radio efectivo** (galaxias con n < 2): dentro del R_e de GalfitM cae el 46 %
  de la luz de la imagen (esperado: 50 %, menos por la PSF); dentro del R_e de
  LSST, solo el 26 %. El radio de media luz medido en la imagen es 1.07× el de
  GalfitM y 1.66× el de LSST. Antes de publicarlo, conviene verificar la
  definición de `sersic_reff_major` en MultiProFit.

---

## Referencias

- Conselice, C. J. (2003), ApJS, 147, 1 — CAS
- Häußler, B. et al. (2013), MNRAS, 430, 330 — GalfitM / MegaMorph
- Hoyos, C. et al. (2011), MNRAS, 411, 2439 — *Residual flux fraction*
- Ivezić, Ž. et al. (2019), ApJ, 873, 111 — LSST
- Lotz, J. M., Primack, J. & Madau, P. (2004), AJ, 128, 163 — Gini–M20
- Rodriguez-Gomez, V. et al. (2019), MNRAS, 483, 4140 — statmorph
- Vika, M. et al. (2014), MNRAS, 444, 3603 — bulbo + disco multibanda con GalfitM

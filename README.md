# Pipeline GalfitM + statmorph + bulbo/disco — LSST DP2 (z < 0.1)

Desde la descarga de las imágenes en el RSP (paso 0) hasta la medición no
paramétrica con **statmorph**, la descomposición **bulbo + disco** con GalfitM
y la **clasificación** de esas descomposiciones, para galaxias con redshift
espectroscópico a z < 0.1.

El pipeline está pensado para trabajar **por partes**: las galaxias se descargan
por lotes, se agrupan en carpetas de trabajo (`galaxias_1_9`, `galaxias_10_22`...)
y cada vez que se suma una carpeta, un solo comando (`actualizar_todo.py`) corre
**solo lo que falta** (pasos 3 a 8). La clasificación bulbo + disco (pasos 9 a 13)
se corre después, a mano, desde `analisis_conjunto/`.

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
| 9 | `analisis_conjunto/` | `leer_output_bd.py` | Lee las salidas bulbo + disco → `galfitm_bd_resultados.csv` |
| 10 | `analisis_conjunto/` | `filtro_bd.py` | Clasifica cada descomposición (`bd_ok`, `bulbo_no_resuelto`, `invertida`...) |
| 11 | `analisis_conjunto/` | `generar_bd_refit.py` + `correr_bd.sh` | (opcional) 2.º ajuste `_bd2` para rescatar invertidas, fallidas y bulbos absurdos |
| 12 | `analisis_conjunto/` | `filtro_una_componente.py` | Reasigna las que en realidad son de una componente (elípticas, discos duplicados) |
| 13 | `analisis_conjunto/` | `inspeccion_visual_bd.py` | Inspección visual interactiva y columna `visual_ok` |

Los pasos 3 a 8 los coordina **`actualizar_todo.py`** (ver [Uso diario](#uso-diario-actualizar_todopy)).
Los pasos 9 a 13 se corren en ese orden (ver [Clasificación bulbo + disco](#clasificación-bulbo--disco-pasos-9-a-13)).

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
    ├── generar_bd_refit.py
    ├── correr_bd.sh
    ├── leer_output_bd.py
    ├── filtro_bd.py
    ├── filtro_una_componente.py
    ├── inspeccion_visual_bd.py
    ├── correr_galfitm_paralelo.sh
    ├── leer_output_galfitm.py
    ├── galfitm-1.4.4-linux-x86_64
    └── analisis_galfitm/                # resultados (tablas y figuras)
        ├── muestra/                     # paso 6
        ├── statmorph_corridas/          # paso 7
        ├── statmorph_resultados.csv     # paso 7
        └── bd/
            ├── galfitm_bd_resultados.csv    # paso 9 (1er ajuste)
            ├── galfitm_bd2_resultados.csv   # paso 11 (2do ajuste, si existe)
            ├── imagenes/                    # figuras de leer_output_bd.py
            └── filtro/                      # pasos 10, 12 y 13
```

Contenido de cada carpeta de galaxia (`galfitm_inputs/<objectId>/`):

| Archivo | Origen |
|---|---|
| `<ID>_sci_<b>.fits`, `<ID>_sig_<b>.fits`, `<ID>_psf_<b>.fits`, `<ID>_mask.fits` | paso 0 |
| `<ID>.feedme` | paso 0 (single-Sérsic) |
| `<ID>_galfitm_out.fits`, `<ID>_galfit.01`, `galfitm_<ID>.log` | paso 3 |
| `<ID>_bd.feedme`, `<ID>_bd.constraints` | paso 8 |
| `<ID>_galfitm_bd_out.fits`, `<ID>_bd_galfit.01`, `galfitm_bd_<ID>.log` | paso 8 |
| `<ID>_bd2.feedme`, `<ID>_bd2.constraints`, `<ID>_galfitm_bd2_out.fits`, `galfitm_bd2_<ID>.log` | paso 11 (solo las reajustadas) |

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
`fallidas_bd.txt`. Con `SUFIJO=_bd2` corre el 2.º ajuste (paso 11).

---

## Clasificación bulbo + disco (pasos 9 a 13)

Todo desde `analisis_conjunto/`, en este orden:

```bash
python3 leer_output_bd.py                     # 9.  lee las descomposiciones
python3 filtro_bd.py                          # 10. clasifica
# 11. (opcional) 2do ajuste de las que fallaron, y volver a leer + filtrar:
#     python3 generar_bd_refit.py --prueba 20 && SUFIJO=_bd2 ./correr_bd.sh 3
#     python3 leer_output_bd.py --sufijo _bd2 && python3 filtro_bd.py
python3 filtro_una_componente.py --solo-ver   # 12. cuenta cuántas reasignaría
python3 filtro_una_componente.py              #     aplica
python3 inspeccion_visual_bd.py               # 13. inspección visual
python3 inspeccion_visual_bd.py --unir        #     agrega visual_ok a las tablas
```

> **Si vuelves a correr `filtro_bd.py`, repite después `filtro_una_componente.py`
> e `inspeccion_visual_bd.py --unir`.** La inspección visual vive en su propio
> archivo, así que no se pierde.

## Paso 9 — Leer salidas bulbo + disco (`leer_output_bd.py`)

Lee el header de cada `<ID>_galfitm_bd_out.fits` de todas las carpetas
`galaxias_*` y crea `analisis_galfitm/bd/galfitm_bd_resultados.csv` (una fila por
galaxia):

- por banda: `B_MAG, B_RE, B_N, B_AR, B_PA` (bulbo), `D_MAG, D_RS, D_RE, D_AR, D_PA`
  (disco, con R_e = 1.678 R_s), `BT`, `MAG_TOT`, `RATIO_RE` (R_e,b / R_e,d);
- `chi2nu`, `params_flag` (parámetros marcados con `*`), `bandas`, `carpeta`.

Si el disco es un Sérsic con n = 1 (2.º ajuste), su `RE` se convierte a R_s.

```bash
python3 leer_output_bd.py                    # solo lee las nuevas
python3 leer_output_bd.py --imagenes 50      # + figuras de 50 galaxias nuevas
python3 leer_output_bd.py --ids analisis_galfitm/bd/filtro/revisar_invertida.csv
                                             # figuras de esa lista (imagenes/revisar_invertida/)
python3 leer_output_bd.py --sufijo _bd2      # 2.º ajuste → galfitm_bd2_resultados.csv
```

## Paso 10 — Clasificación (`filtro_bd.py`)

Filtro lógico basado en **Vika et al. (2014)** y adaptado a la resolución de LSST,
pensado para la relación masa–tamaño de bulbos y discos. Usa la banda **r** como
referencia (i si falta): R_e, n, b/a y PA son constantes en λ, así que basta una.

Cada galaxia recibe **una sola categoría**, la del primer criterio que cumple, en
este orden:

| # | Categoría | Condición | Qué se usa |
|---|---|---|---|
| 1 | `fallida` | R_e, n, b/a del bulbo o R_s, b/a del disco marcados con `*`; χ²/ν fuera de 0.5–5; \|m(B+D) − m(single-Sérsic)\| > 0.3 mag; B/T vacío; o (2.º ajuste) solución pegada al borde de las restricciones | descartar o reajustar |
| 2 | `perturbada` | asimetría de statmorph A > 0.35 (Conselice 2003) | descartar (fusión, interacción, vecino) |
| 3 | `invertida` | R_e,b ≥ R_e,d: las componentes intercambiaron papeles | reajustar (paso 11) |
| 4 | `esferoide_puro` | B/T ≥ 0.94 | single-Sérsic |
| 5 | `disco_puro` | B/T ≤ 0.06 | single-Sérsic como disco |
| 6 | `disco_absurdo` | R_e,d > 3 × R_e del single-Sérsic (el disco ajusta cielo o halo) | descartar |
| 7 | `posible_barra` | n_b < 0.7, (b/a)_b < 0.5, \|ΔPA\| > 20° y R_e,b < R_e,d | excluir (comparación con muestras sin barras) |
| 8 | `bulbo_aguja` | (b/a)_b < 0.2: barra, disco interno o artefacto | excluir |
| 9 | `disco_de_canto` | (b/a)_d < 0.3: B+D no confiable de canto | excluir |
| 10 | `bulbo_absurdo` | n_b < 0.3 o n_b > 8 | reajustar (paso 11) |
| 11 | `bulbo_no_resuelto` | R_e,b < 1 FWHM de la PSF | B/T y disco sí; tamaño del bulbo como **límite superior** |
| 12 | `bd_ok` | pasa todo | **muestra para el análisis** |

Detalles:

- **FWHM de la PSF**: se mide en la PSF de cada galaxia (`<ID>_psf_<b>.fits`) como
  2√(A/π), con A el área sobre la mitad del máximo. Se guarda en `psf_fwhm.csv` y
  solo se miden las nuevas.
- **Asimetría**: se toma de `statmorph_resultados.csv` en r, i o g (la primera con
  `flag ≤ 1` y S/N ≥ 2.5). Si no existe la tabla, el corte no se aplica.
- Cada criterio queda además en una columna booleana `corte_*`, para poder
  recombinarlos sin volver a correr nada.
- Si hay z, calcula `RE_B_KPC`, `RE_D_KPC` y `FWHM_KPC`.
- Si existe `galfitm_bd2_resultados.csv`, clasifica los dos ajustes y se queda,
  para cada galaxia, con el de **mejor categoría** (en empate, el 1.º). La columna
  `ajuste` (`bd` / `bd2`) dice de cuál viene, y en pantalla sale una tabla de
  rescate (categoría con el 1.º → con el 2.º).

Variantes (cada una en su subcarpeta, p. ej. `filtro/ba_min_0.5/`; `psf_fwhm.csv`
se comparte):

```bash
python3 filtro_bd.py --fwhm-min 1.5     # resolución más exigente
python3 filtro_bd.py --ba-min 0.5       # además disco poco inclinado → si no, 'disco_inclinado'
python3 filtro_bd.py --banda i          # otra banda de referencia
```

Salidas en `analisis_galfitm/bd/filtro/`:

| Archivo | Contenido |
|---|---|
| `galfitm_bd_filtrado.csv` | todas, con `categoria_bd`, `ajuste` y columnas `corte_*` |
| `muestra_bd_ok.csv` | solo `bd_ok` |
| `revisar_<categoria>.csv` | 20 al azar de cada categoría (para `--ids`) |
| `psf_fwhm.csv` | caché de la FWHM de la PSF |
| `22_filtro_bd.png` | diagnósticos: B/T, R_e,b/R_e,d, n_b, R_e,b/FWHM, R_e,b vs R_e,d, conteo por categoría |

En pantalla: conteo por categoría, cuántas fallan cada corte, medianas de `bd_ok`,
cuántas `bd_ok` tienen R_e,b ≥ 1.5 y ≥ 2 FWHM, y la tabla categoría × tipo
morfológico (ETG/LTG/TG/OG).

## Paso 11 — 2.º ajuste (`generar_bd_refit.py`, opcional)

Rehace el ajuste de las `invertida`, `bulbo_absurdo` y `fallida` (con
`--color-malo` también las que tuvieron masa con el color de la galaxia o masa
inconsistente). Usa archivos propios (`_bd2`), así que **el 1.º ajuste no se toca**.

| Problema del 1.º ajuste | Qué cambia en el 2.º |
|---|---|
| Invertidas | Restricción de cociente `1/2 re 0.02 0.9` (R_e,b ≤ 0.9 R_e,d). Para usarla, el disco pasa a Sérsic con n = 1 fijo (mismo perfil que el expdisk, pero con `re`) |
| n del bulbo absurdo | n_b en 0.3–8, empieza en n = 2 |
| Componentes que desaparecen en una banda | magnitudes como polinomio suave en λ (lineal con ≤ 4 bandas, cuadrático con ≥ 5) y `1-2 mag -3 3` (0.06 < B/T < 0.94 en cada banda) |
| Bulbos de menos de un píxel | R_e,b ≥ 0.5 px |
| Inicio que favorece la inversión | bulbo chico y débil (0.2 R_e,ss, m_ss + 1.5), disco grande y brillante (1.2 R_e,ss, m_ss + 0.4) |

```bash
python3 generar_bd_refit.py --prueba 20     # primero 20 al azar
SUFIJO=_bd2 ./correr_bd.sh 3
python3 leer_output_bd.py --sufijo _bd2 --ver 20
python3 filtro_bd.py
python3 generar_bd_refit.py                 # si la prueba salió bien: todas
```

En la prueba, revisar que el log `galfitm_bd2_<ID>.log` no tenga errores al leer
las restricciones y que en los resultados R_e,b/R_e,d ≤ 0.9 y B/T por banda quede
en 0.06–0.94. Una solución pegada a esos bordes (cociente > 0.88, n en el límite
o |Δm| > 2.95) la forzó la restricción y no los datos: `filtro_bd.py` la marca
como `fallida`.

## Paso 12 — Galaxias de una componente (`filtro_una_componente.py`)

Revisa las `bd_ok` y `bulbo_no_resuelto` y reasigna las que en realidad se
describen con **una sola componente**, usando categorías que el resto del
pipeline ya entiende:

| Criterio | Condición | Pasa a |
|---|---|---|
| `sin_cruce` | el bulbo domina el brillo superficial en todos los radios entre 0.5 FWHM y 3 R_e,d (perfiles 1D de cada componente): el "disco" solo corrige el perfil de una elíptica | `esferoide_puro` |
| `disco_duplicado` | n_b < 1.5, R_e,b/R_e,d > 0.6, \|Δ(b/a)\| < 0.1 y \|ΔPA\| < 10° (el PA no se exige si las dos tienen b/a > 0.85): el mismo disco repartido en dos exponenciales | `disco_puro` |

- Agrega a `galfitm_bd_filtrado.csv`: `categoria_bd_antes`, `una_comp_motivo`,
  `corte_sin_cruce`, `corte_disco_duplicado` y `r_disco_domina_arcsec`.
- Reescribe `muestra_bd_ok.csv` y deja `revisar_sin_cruce.csv` y
  `revisar_disco_duplicado.csv` (20 al azar) para verificarlas a ojo.
- Se puede correr varias veces: siempre parte de lo que dejó `filtro_bd.py`.
- Umbrales al principio del script: `R_MAX_RED`, `DUP_N_MAX`, `DUP_RATIO_MIN`,
  `DUP_DQ_MAX`, `DUP_DPA_MAX`, `DUP_Q_REDONDA`.

## Paso 13 — Inspección visual (`inspeccion_visual_bd.py`)

Ventana interactiva para aceptar o rechazar cada descomposición. Por defecto
muestra `bd_ok`, `bulbo_no_resuelto` y `disco_puro` (sin las reasignadas
automáticamente en el paso 12), ordenadas por categoría.

Cada figura muestra: imagen, modelo B+D, residuo y residuo relativo; bulbo,
disco y perfil radial en anillos elípticos (dato, modelo, bulbo y disco);
elipses del R_e del bulbo (rojo) y del disco (azul), con la misma escala de
brillo en imagen, modelo y componentes.

| Tecla | Acción |
|---|---|
| `a` / `1` | aceptar |
| `x` / `0` | rechazar |
| `d` / `2` | dudosa |
| `b` `o` `v` `s` `e` `i` | etiquetas: barra, anillo, vecino, estructura, artefacto, componentes intercambiadas (antes de decidir) |
| `→` / `←` | moverse sin decidir |
| `n` | próxima sin revisar |
| `f` | cambiar de banda |
| `+` / `-` | acercar / alejar |
| `q` | salir |

El progreso se guarda después de cada decisión en
`analisis_galfitm/bd/filtro/inspeccion_visual.csv` (escritura atómica): al volver
a abrirlo empieza en la primera sin revisar.

```bash
python3 inspeccion_visual_bd.py --revisar dudosa     # volver a ver las dudosas
python3 inspeccion_visual_bd.py --orden azar
python3 inspeccion_visual_bd.py --ids analisis_galfitm/bd/filtro/revisar_sin_cruce.csv
python3 inspeccion_visual_bd.py --png 50             # sin ventana: guarda 50 PNG
python3 inspeccion_visual_bd.py --resumen            # cuántas llevas, por decisión y etiqueta
python3 inspeccion_visual_bd.py --unir               # agrega las columnas a las tablas
```

`--unir` agrega `visual`, `visual_ok` (1 = aceptada, 0 = rechazada, 2 = dudosa,
vacío = sin revisar) y `visual_etiquetas` a `galfitm_bd_filtrado.csv` y a
`componentes_masas.csv`, guardando la primera vez una copia `.bak` de cada tabla.

---

## Uso diario: `actualizar_todo.py`

Desde `analisis_conjunto/`. Recorre las carpetas `galaxias_*` (junto a esta
carpeta o un nivel arriba) y corre solo lo pendiente:

1. single-Sérsic de galaxias sin salida (no reintenta las fallidas) y
   `leer_output_galfitm.py` si hay salidas nuevas;
2. `combinar_carpetas.py` y `procesar_muestra.py` (cachés: solo mide lo nuevo);
3. statmorph y bulbo + disco **solo** para las A + B que no los tienen;
4. une los resultados de statmorph.

La clasificación bulbo + disco (pasos 9 a 13) **no** la corre: después de
`actualizar_todo.py`, correrla a mano (ver arriba).

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

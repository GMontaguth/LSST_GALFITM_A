# Pipeline GalfitM + statmorph + bulbo/disco — LSST DP2 (z < 0.1)

Desde la descarga de las imágenes en el RSP (paso 0) hasta la medición no
paramétrica con **statmorph**, la descomposición **bulbo + disco** con GalfitM
y la **clasificación** de esas descomposiciones, para galaxias con redshift
espectroscópico a z < 0.1.

El pipeline está pensado para trabajar **por partes**: las galaxias se descargan
por lotes, se agrupan en carpetas de trabajo (`galaxias_1_9`, `galaxias_10_22`...)
y cada vez que se suma una carpeta, un solo comando (`actualizar_todo.py`) corre
**solo lo que falta** (pasos 3 a 8). La clasificación bulbo + disco (pasos 9 a 14)
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
| 14 | `analisis_conjunto/` | Claude Code (terminal) | Clasificación visual asistida: calibración a ciegas contra la humana y clasificación del resto |

Los pasos 3 a 8 los coordina **`actualizar_todo.py`** (ver [Uso diario](#uso-diario-actualizar_todopy)).
Los pasos 9 a 14 se corren en ese orden (ver [Clasificación bulbo + disco](#clasificación-bulbo--disco-pasos-9-a-14)).

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
            └── filtro/                      # pasos 10, 12, 13 y 14
                ├── inspeccion_visual.csv         # decisiones humanas
                ├── inspeccion_visual_claude.csv  # decisiones de Claude
                └── inspeccion_png/               # figuras (--png)
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

## Clasificación bulbo + disco (pasos 9 a 14)

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
# 14. clasificación asistida con Claude Code (ver abajo), y después:
python3 inspeccion_visual_bd.py --unir --final   # visual_ok_final (humana > Claude > auto)
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
muestra `bd_ok`, `bulbo_no_resuelto` y `disco_puro`, ordenadas por categoría.
Omite las que `filtro_una_componente.py` ya reasignó como una componente (avisa
cuántas son); con `--incluir-auto` también las muestra.

Cada figura muestra: imagen, modelo B+D, residuo y residuo relativo; bulbo,
disco y perfil radial en anillos elípticos (dato, modelo, bulbo y disco);
elipses del R_e del bulbo (rojo) y del disco (azul), con la misma escala de
brillo en imagen, modelo y componentes.

| Tecla | Acción |
|---|---|
| `a` / `1` | aceptar |
| `x` / `0` | rechazar |
| `d` / `2` | dudosa |
| `u` / `3` | una componente (el ajuste no es malo, pero no hacen falta dos componentes) |
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
3 = una componente, vacío = sin revisar) y `visual_etiquetas` a
`galfitm_bd_filtrado.csv` y a `componentes_masas.csv`, guardando la primera vez
una copia `.bak` de cada tabla.

## Paso 14 — Clasificación visual asistida con Claude (Claude Code)

Para no revisar a mano las ~1500 descomposiciones, se usa **Claude Code** en la
terminal como segundo clasificador. Trabaja directamente en `analisis_conjunto/`:
genera las figuras con el inspector, las mira una por una y escribe sus
decisiones con los mismos códigos y etiquetas. **La clasificación humana manda:**
Claude escribe en un archivo aparte y nunca modifica `inspeccion_visual.csv`.

```bash
curl -fsSL https://claude.ai/install.sh | bash      # instalar (una vez)
cd "/media/gissel/Sem Título/DP2_LSST/Galfitm/analisis_conjunto"
claude                                              # iniciar sesión y aceptar la carpeta
```

### Procedimiento

1. **Criterios comunes.** Claude recibe los mismos criterios que la inspección
   humana (abajo) y lee el código de `inspeccion_visual_bd.py` para entender cada
   panel de la figura.
2. **Calibración a ciegas.** Clasifica las ~100 galaxias que ya revisé a mano
   **sin leer** `visual`, `visual_ok` ni `visual_etiquetas` (usa el CSV solo para
   sacar la lista de `objectId`). Las figuras se generan con `--ids` y
   `--incluir-auto`, para que no falte ninguna. Después compara: matriz de
   confusión, **κ de Cohen** y lista de desacuerdos con su razonamiento. También
   informa cuántas de mis galaxias reasignó automáticamente
   `filtro_una_componente.py` y si mi decisión coincide con la automática.
   **Se detiene aquí** hasta que yo apruebe.
3. **Resto de la muestra.** Si la calibración es buena, clasifica las demás (sin
   `--incluir-auto`: las reasignadas automáticamente conservan su clase) y guarda
   en **`analisis_galfitm/bd/filtro/inspeccion_visual_claude.csv`**, con las
   mismas columnas más `clasificador = claude` y una nota breve por galaxia.
   Guarda por lotes, así que se puede retomar si se corta.
4. **Control de lo automático.** Revisa 20 al azar de las reasignadas por
   `filtro_una_componente.py` (con `--incluir-auto`) e informa si está de acuerdo
   con la clasificación de una componente.
5. **Revisión humana.** Reviso las dudosas de Claude y una fracción al azar de
   sus decisiones.

### Criterios (los mismos para la inspección humana y la de Claude)

| `visual_ok` | Decisión | Criterio |
|---|---|---|
| 1 | aceptada | bulbo y disco son componentes reales: en el perfil radial el disco domina en las afueras (las curvas se cruzan), el residuo no tiene estructura fuerte, y R_e y b/a tienen sentido físico |
| 3 | una componente | el ajuste no es malo, pero no hacen falta dos componentes: el bulbo domina en todos los radios (elípticas, componentes redondas y parecidas) o las dos describen el mismo disco (n_b bajo, R_e, b/a y PA parecidos) |
| 0 | rechazada | ajuste malo: barra o anillo fuerte mal absorbido, vecino brillante que contamina, o componentes intercambiadas (el "bulbo" más grande y extendido que el disco) |
| 2 | dudosa | la figura no permite decidir |

- **No** se rechaza solo por el patrón rojo/azul en el centro (< 2″, típico de la
  PSF), por manchas de vecinos enmascarados ni por el ruido en el borde.
- Etiquetas (se combinan con cualquier decisión): `b` barra, `o` anillo,
  `v` vecino, `s` estructura, `e` artefacto, `i` componentes intercambiadas.

### Limitaciones

- Claude clasifica a partir de la **figura PNG**, no del FITS: lo que no se vea
  con ese estiramiento de brillo tampoco lo ve.
- No es un clasificador validado. En los casos claros coincide casi siempre; en
  los límites (barras débiles, anillos sutiles) puede diferir, como dos
  inspectores humanos. Por eso la calibración no es opcional.
- Si una decisión humana y una de Claude no coinciden, gana la humana.

### Resultados de la calibración (1.ª ronda, a ciegas)

101 galaxias revisadas a mano; Claude las clasificó con figuras generadas con
`--limpia` (sin la decisión humana impresa) y solo después leyó las columnas
humanas. Detalle por galaxia en `analisis_galfitm/bd/filtro/calibracion_claude.csv`.

| humana ↓ / Claude → | 1 aceptada | 3 una comp. | 0 rechazada | 2 dudosa | total |
|---|---|---|---|---|---|
| 1 aceptada | 36 | 8 | 9 | 2 | 55 |
| 3 una componente | 6 | 31 | 0 | 0 | 37 |
| 0 rechazada | 1 | 0 | 4 | 0 | 5 |
| 2 dudosa | 3 | 0 | 1 | 0 | 4 |
| total | 46 | 39 | 14 | 2 | 101 |

| Comparación | Acuerdo | κ de Cohen |
|---|---|---|
| Las cuatro clases | 71/101 (70 %) | 0.51 |
| Sin las dudosas de ninguno de los dos | 71/95 | 0.56 |
| Aceptada vs. el resto | 72/101 | 0.43 |

Desacuerdos (30), en tres grupos sistemáticos:

1. **Barras (9):** aceptadas a mano con etiqueta `b`, rechazadas por Claude
   porque el "bulbo" es la propia barra o deja un residuo estructurado.
2. **Bulbos de n bajo con la forma del disco (8):** aceptadas a mano, una
   componente para Claude.
3. **Bulbos compactos y débiles en espirales (6):** una componente a mano,
   aceptadas por Claude.

Solo 1 de las 101 había sido reasignada por `filtro_una_componente.py`
(`sin_cruce`), y las tres clasificaciones coinciden. En total hay 139
reasignadas automáticamente (75 `disco_duplicado` + 64 `sin_cruce`).

Claude marca más etiquetas que la inspección humana (estructura 35 vs. 16, barra
20 vs. 14, incluidas las 14 humanas).

> El κ reportable es el de esta ronda a ciegas (0.51). Si se cambian las reglas,
> la concordancia recalculada sobre estas mismas 101 (p. ej. κ ≈ 0.63 aceptando
> las barras) ya no es ciega: hay que hacer una 2.ª ronda con galaxias nuevas.

Cambio en el código durante la calibración: `inspeccion_visual_bd.py` fallaba
con galaxias cuyo borde es todo ceros (error en las isofotas); se agregó un
respaldo de tres líneas para estimar el ruido.

### Mensajes a Claude Code

Registro textual de las instrucciones, en el orden en que se enviaron, para que
el procedimiento sea reproducible (y para el apéndice del paper).

**Mensaje 1 — instrucciones y calibración a ciegas** (enviado el 2026-10-13).

> Nota: el mensaje dice 75 galaxias reasignadas automáticamente; en la tabla
> eran 139 (75 `disco_duplicado` + 64 `sin_cruce`), como informó Claude Code.

```
You will act as my visual classifier for GalfitM bulge+disc decompositions.

Context:
- The inspection tool is inspeccion_visual_bd.py. Read its code to understand what each
  figure panel shows (image, model, residual, residual/model, bulge, disc and radial
  profile) and its options (--png, --ids, --incluir-auto, etc.).
- My decisions are stored in analisis_galfitm/bd/filtro/inspeccion_visual.csv (I have
  inspected about 100 galaxies). During calibration, use this file ONLY to get the list
  of objectIds I inspected; do NOT read the visual, visual_ok or visual_etiquetas
  columns until you have finished classifying them.
- filtro_una_componente.py has already classified 75 galaxies as single component
  automatically; the inspector skips them by default (--incluir-auto shows them).
- You do not need any notebook.

Classification codes (visual_ok):
- 1 accepted: the bulge and disc are real components. In the radial profile the disc
  dominates in the outskirts (the profiles cross), the residual shows no strong
  structure, and the effective radii and axis ratios are physically meaningful.
- 3 single component: the fit is not bad, but the galaxy does not require two
  components. Examples: the bulge dominates at all radii (typical of ellipticals, with
  both components round and similar), or both components describe the same disc (low
  n_b, similar R_e, b/a and PA).
- 0 rejected: the fit is poor. Examples: a strong bar or ring that is poorly absorbed
  by the model, a bright neighbour contaminating the fit, or swapped components (the
  "bulge" is larger and more extended than the disc).
- 2 uncertain: the figure does not allow a decision.
- Do not reject a galaxy based only on the red/blue pattern in the centre (< 2″,
  typical of PSF mismatch), on residuals from masked neighbours, or on noise at the
  outer edge.
- Flags (can be combined with any decision): b bar, o ring, v neighbour, s structure,
  e artefact, i swapped components.

Plan:
1. CALIBRATION: generate figures for ALL the objectIds I inspected, using --ids with my
   list and --incluir-auto so none of them is left out. Inspect each figure yourself and
   classify them blind to my decisions. Then compare with mine: confusion matrix,
   Cohen's kappa, and the list of disagreements with your reasoning. Also report how
   many of my inspected galaxies were auto-classified as single component and whether
   my decision agrees with the automatic one. STOP and show me this before continuing.
2. If I approve, generate figures for the remaining galaxies (without --incluir-auto),
   classify them, and save them to a SEPARATE file,
   analisis_galfitm/bd/filtro/inspeccion_visual_claude.csv, with the same columns plus
   clasificador=claude and a short note per galaxy. Do NOT modify my
   inspeccion_visual.csv.
3. The auto-classified galaxies keep their automatic class, but inspect a random
   subsample of 20 of them (with --incluir-auto) and report whether you agree with the
   automatic single-component classification.
4. Save progress in batches so the work can be resumed if interrupted.
```

**Mensaje 2 — primeras reglas tras la calibración.** Resuelve los tres grupos
de desacuerdo de la 1.ª ronda y pide una 2.ª ronda a ciegas antes de clasificar
el resto.

```
Rule changes after calibration:
1. Bars: accept with flag b if the bulge is round/aligned and the bar is a separate
   residual feature. Reject if the "bulge" component itself is the bar (elongated,
   b/a < ~0.4, misaligned with the disc).
2. Compact weak bulges: accept only if the bulge dominates out to a radius larger
   than 1 PSF FWHM (use r_disco_domina_arcsec from galfitm_bd_filtrado.csv). If the
   profiles cross inside 1 FWHM, classify as single component.
3. Low-n disc-shaped bulge: accept if R_e,b/R_e,d < 0.5 (possible pseudobulge).
   Single component only if R_e are similar AND the "bulge" holds > 50% of the light.
Do NOT start step 2 yet. First I will classify 50 new random galaxies myself; then
classify those same 50 blind with these rules and report the confusion matrix and
kappa again.
```

**Mensaje 3 — corrección del filtro y revisión de las reglas.** Claude Code
encontró que `filtro_una_componente.py` tomaba `RE_B_FWHM` como FWHM de la PSF
(ver *Correcciones de código*, abajo).

```
Thanks, good catches. Before the second round:
1. Fix the FWHM column bug in filtro_una_componente.py so it uses FWHM_PSF, rerun it,
   and report what changed (sin_cruce, disco_duplicado, r_disco_domina_arcsec).
   Do not decide those 12 galaxies by eye; use the corrected radius.
2. Revise rule 2: single component only if the profiles cross inside 1 FWHM AND
   B/T < 0.2. If B/T >= 0.2 and they cross inside 1 FWHM, accept it as an unresolved
   bulge (B/T usable, bulge size is an upper limit).
3. Rule 1 (bars) is a visual criterion, not a threshold: reject only when the figure
   shows that the "bulge" component is the bar itself.
4. disco_puro galaxies are single component by definition: for them, only decide
   whether the single-component fit is good (3) or bad (0). Rule 2 does not apply.
5. List the calibration galaxies where the revised rules disagree with my decision,
   so I can review them and either change my decision or the rule.
The second round will use the default mixed selection. Still do not start step 2.
```

**Mensaje 4 — decisiones sobre las reglas.** Antes se guardó una copia de las
decisiones originales (`inspeccion_visual_ronda1.csv`) y se revisaron a mano las
25 galaxias en conflicto.

```
Decisions on the rules:
1. Rule 3: keep R_e,b/R_e,d >= 0.5 as "similar R_e". Low-n bulges (n_b < ~0.6) with
   similar R_e and B/T > 0.5 are single component.
2. Rule 2: the "unresolved bulge" branch applies only to compact bulges
   (R_e,b/R_e,d < 0.5). Extended "bulges" crossing inside 1 FWHM: decide visually,
   usually as a split disc (single component).
3. Rule 2 takes precedence over rule 3.
4. Weak needle components in discs with B/T < ~0.11: single component, not rejected.
I reviewed the 25 calibration galaxies and updated my decisions (my original round-1
decisions are saved in inspeccion_visual_ronda1.csv; use that file for the round-1
kappa).
```

**Mensaje 5 — barras, reglas 3 y 4 como guía, y arreglo del inspector.**

```
Correction about bars: I do reject for bars sometimes (739160573106737611 and
739163802922149219). My rule: if the "bulge" component is the bar, accept with flag b
(these galaxies will be excluded from the bulge analysis). If the "disc" component
fits the bar, or the residual shows the fit is broken, reject.
I kept 739156656096571847 and 772395442360045399 as accepted, so do not apply rules
3 and 4 as hard thresholds; use them as guidance and decide those cases from the figure.
Also patch inspeccion_visual_bd.py to re-read and merge before each save (by objectId,
keeping the most recent decision per galaxy), tested on a copy of the file.
```

**Mensaje 6 — 2.ª ronda a ciegas.**

```
I have finished classifying 50 new galaxies. Run the blind second round:
1. Take as the round-2 set the objectIds in inspeccion_visual.csv that are not in the
   calibration list. Tell me how many you find and their categoria_bd.
2. Generate their figures with --limpia so my decisions are not printed on them.
   Do NOT read my visual, visual_ok or visual_etiquetas columns for these galaxies
   until you have classified all of them.
3. Classify them with the final rules.
4. Then compare: confusion matrix, Cohen's kappa, agreement excluding uncertain, and
   accepted vs. rest, plus the list of disagreements with your reasoning. Report it
   also split by categoria_bd (bd_ok, bulbo_no_resuelto, disco_puro).
Save your round-2 calls to calibracion_claude_ronda2.csv. Still do not start step 2.
```

**Mensaje 7 — aprobación y paso 2.**

```
Round 2 approved. Before step 2:
- 780665900104701331 and 792329931769003310 were slips; I will correct them to single
  component (3). Keep the rules as written.
- Rule 2 tolerance: if the crossing radius is within ±10% of 1 FWHM, decide from the
  figure instead of the threshold.
- Reject only when the fit is plainly unusable. Otherwise keep your decision and use
  flags (v neighbour, e artefact, s structure).
- Flag b only when the bar is evident: galaxies with flag b will be excluded from the
  bulge analysis, so over-flagging removes valid galaxies.
Start step 2 on the remaining galaxies, in saved batches of ~50, writing only to
inspeccion_visual_claude.csv (clasificador=claude, short note per galaxy). Do not
modify my inspeccion_visual.csv. At the end, give me:
- all your rejected and uncertain galaxies,
- all galaxies with flag b,
- a random sample of 50 of the rest (fixed seed),
as separate CSVs I can open with inspeccion_visual_bd.py --ids, plus the auto-classified
check of 20 (step 3 of the original plan).
```

**Mensaje 8 — 3.ª ronda y columna final.** (Primero se revisó a mano la muestra
de 50 al azar, solo con el inspector, sin abrir el CSV con las decisiones de Claude.)

```
I have reviewed your step-2 samples in the inspector; my decisions are now in
inspeccion_visual.csv. Some lists I may have reviewed only partly, so first tell me
how many galaxies of each review file are now in my CSV.
1. Round 3: compare my decisions with yours on revisar_claude_muestra50.csv: confusion
   matrix, Cohen's kappa (all classes and accepted vs. rest), disagreements with your
   reasoning, and agreement split by the six parallel passes.
2. For the bar-flagged galaxies and the 13 rejected/uncertain that I reviewed: how many
   I changed, and in which direction. For bars, also whether I kept or removed flag b.
3. Extend inspeccion_visual_bd.py --unir (keep the current behaviour as default) with
   an option that writes:
   - visual_ok_final: my decision if it exists, otherwise Claude's, otherwise the
     automatic single-component class from filtro_una_componente.py (3);
   - visual_fuente: 'humana', 'claude' or 'auto';
   - visual_etiquetas_final with the same priority.
   Test it on a copy of the tables first.
Do not modify inspeccion_visual.csv.
```

**Mensaje 9 — disco_puro, origen de las "una componente" y lista de barras.**

```
I reviewed the 13 rejected/uncertain galaxies (revisar_claude_rechazadas_dudosas.csv);
my decisions are in inspeccion_visual.csv. Before merging:
1. disco_puro: when I accept a disco_puro, I mean "the single-component fit is good".
   In --final, map visual_ok = 1 to 3 for disco_puro galaxies (mine and yours).
2. Rule 3 split-disc branch: across rounds I consistently accept low-n bulges with
   similar R_e and B/T > 0.5, and you call them single. Count how many of your 1,331
   calls are single component mainly because of this branch (n_b < ~1,
   R_e,b/R_e,d >= 0.5, B/T > 0.5), and how many are single by rule 2 with
   B/T < 0.2. Give me both lists as CSVs, plus 10 random IDs from the first list
   that I can review. Do not reclassify anything yet.
3. Did you implement the re-read-and-merge-before-save patch in
   inspeccion_visual_bd.py? Confirm.
4. Make revisar_claude_barra20.csv with 20 random galaxies (fixed seed) from
   revisar_claude_barra.csv, so I can review them.
5. Compare my decisions on the 13 rejected/uncertain with yours: how many I changed
   and in which direction.
Do not run --unir on the real tables yet.
```

**Mensaje 10 — disco partido según el criterio humano y muestra "tipo elíptica".**

```
Thanks.
1. Split-disc branch: follow my criterion. Reclassify those 17 galaxies (those without
   a decision of mine) to accepted in inspeccion_visual_claude.csv, with a note
   "reclassified: split-disc rule follows human criterion". Report how many changed.
2. Elliptical-like group (B/T > 0.5, n_b >= 2, your single calls): make
   revisar_claude_eliptica15.csv with 15 random galaxies (fixed seed) that I have not
   reviewed, so I can check whether some are S0s with a real disc.
Do not run --unir on the real tables yet.
```

**Mensaje 11 — prueba del criterio para el grupo "tipo elíptica".**

```
I have now reviewed revisar_claude_eliptica15.csv (saved in inspeccion_visual.csv).
1. Elliptical-like group: use the disc-excess criterion (accept if the disc is more
   than 0.6 mag above the bulge at 2 R_e,d), not B/T < 0.63. A B/T cut would select on
   a quantity tied to bulge mass and could bias the bulge mass-size slope.
   First test it on my 15 new decisions (they were not used to tune the threshold)
   and report the agreement. If it is reasonable (~80% or more), apply it in both
   directions to your unreviewed calls in this group, with a note on each change,
   and report how many changed. Keep a copy of the file before changing it.
2. Bars: do not change anything. In the analysis I will test the bulge mass-size
   relation with and without flag-b galaxies.
3. Then, if everything is consistent, run inspeccion_visual_bd.py --unir --final on
   the real tables (it makes .bak copies) and report the final counts per
   visual_fuente and visual_ok_final.
```

**Mensaje 12 — unión final.** (El criterio del mensaje 11 no pasó la prueba y no
se aplicó.)

```
Agreed: leave the elliptical-like group unchanged. Run inspeccion_visual_bd.py
--unir --final on the real tables now and report:
- the final counts per visual_fuente and visual_ok_final, for both tables;
- how many galaxies with visual_ok_final = 1 have flag b (for the robustness test);
- confirmation that the .bak copies were created.
```

### Resultados de la 2.ª ronda (a ciegas, reglas finales)

52 galaxias nuevas, clasificadas primero a mano (`--orden azar --semilla 7`) y
después por Claude con figuras `--limpia`: 29 `bd_ok`, 8 `bulbo_no_resuelto`,
15 `disco_puro`, ninguna reasignada automáticamente. Detalle en
`calibracion_claude_ronda2.csv`.

| humana ↓ / Claude → | 1 aceptada | 3 una comp. | 0 rechazada | 2 dudosa | total |
|---|---|---|---|---|---|
| 1 aceptada | 16 | 4 | 1 | 0 | 21 |
| 3 una componente | 2 | 27 | 2 | 0 | 31 |
| 0 rechazada | 0 | 0 | 0 | 0 | 0 |
| 2 dudosa | 0 | 0 | 0 | 0 | 0 |
| total | 18 | 31 | 3 | 0 | 52 |

| Comparación | Acuerdo | κ de Cohen |
|---|---|---|
| Las cuatro clases | 43/52 (83 %) | **0.66** |
| Aceptada vs. el resto | 45/52 | 0.71 |
| Solo `bd_ok` (29) | 23 (79 %) | 0.60 |
| Solo `bulbo_no_resuelto` (8) | 7 (88 %) | — (casi todas en una clase) |
| Solo `disco_puro` (15) | 13 (87 %) | — (casi todas en una clase) |

Evolución del acuerdo:

| Ronda | Galaxias | Reglas | κ | ¿Ciego? |
|---|---|---|---|---|
| 1 | 101 (`bd_ok`) | solo criterios generales | 0.51 | sí |
| 1, recalculada | 101 | reglas finales, decisiones humanas revisadas | 0.71 | **no** |
| 2 | 52 (mezcla real) | reglas finales | **0.66** | sí |

Diferencias que quedaron: Claude rechaza más (3 contra 0, en casos que a mano se
conservan con etiqueta), marca más barras (11 contra 5) y menos estructura
(14 contra 21). Dos decisiones humanas contradecían las reglas y se corrigieron
como descuidos. Esto se ajustó en el mensaje 7.

### Reglas finales (inspección humana y de Claude)

Se suman a la tabla de criterios de arriba.

| Regla | Caso | Decisión |
|---|---|---|
| 1. Barras | el "bulbo" es la barra (alargado, desalineado del disco) | aceptada + `b` |
| | el "disco" ajusta la barra, o el residuo muestra que el ajuste está roto | rechazada |
| 2. Bulbos compactos | perfiles que se cruzan dentro de 1 FWHM (`r_disco_domina_arcsec`), B/T < 0.2 | una componente |
| | ídem con B/T ≥ 0.2 y R_e,b/R_e,d < 0.5 | aceptada como bulbo no resuelto (R_e,b es límite superior) |
| | ídem con R_e,b/R_e,d ≥ 0.5 | se decide con la figura (normalmente disco partido: una componente) |
| | cruce a ±10 % de 1 FWHM | se decide con la figura |
| 3. Bulbos de n bajo con la forma del disco (guía) | R_e,b/R_e,d < 0.5 | aceptada (posible pseudobulbo) |
| | R_e,b/R_e,d ≥ 0.5 y B/T > 0.5 | una componente (mismo disco partido en dos) |
| 4. Componentes aguja débiles (guía) | B/T < ~0.11 en un disco | una componente, no rechazada |
| 5. `disco_puro` | ya es una componente | solo 3 (ajuste bueno) o 0 (malo) |
| Prioridad | | la regla 2 manda sobre la 3 |
| Rechazo | | solo si el ajuste es claramente inutilizable; si no, decisión + etiqueta (`v`, `e`, `s`) |
| Etiqueta `b` | | solo con barra evidente |

**En el análisis**, las galaxias con etiqueta `b` se excluyen de la relación
masa–tamaño de los bulbos (su R_e y n pueden ser los de la barra), lo que además
hace la muestra comparable con Méndez-Abreu et al. (2021), que no tiene barras.

### Paso 2: clasificación de Claude

Claude clasificó las 1381 galaxias restantes (`bd_ok`, `bulbo_no_resuelto` y
`disco_puro`, sin las reasignadas automáticamente ni las ya revisadas a mano) en
seis tandas paralelas con las mismas reglas, guardando por lotes de ~50
(`paso2_claude.py`, `inspeccion_visual_claude.csv`).

| `categoria_bd` | aceptada | una comp. | dudosa | rechazada | total |
|---|---|---|---|---|---|
| `bd_ok` | 347 | 375 | 7 | 1 | 730 |
| `bulbo_no_resuelto` | 296 | 129 | 1 | 1 | 427 |
| `disco_puro` | 0 | 221 | 0 | 3 | 224 |
| total | 643 | 725 | 8 | 5 | 1381 |

- Regla 2: 215 galaxias con cruce dentro de 1 FWHM y B/T < 0.2 → una componente;
  81 bulbos compactos con B/T ≥ 0.2 → aceptadas; 52 extendidos → una componente
  (decididos con la figura); 102 en la franja de ±10 % → decididas con la figura.
- 101 galaxias (7 %) con `b`. Solo 5 rechazos, todos claramente inutilizables.
- Control de las automáticas: 20 al azar de las 139, todas de acuerdo con
  "una componente".

### Control humano del paso 2

**3.ª ronda: 50 al azar de las 1381** (revisadas a mano solo con el inspector,
sin ver las decisiones de Claude):

| humana ↓ / Claude → | 1 aceptada | 3 una comp. | total |
|---|---|---|---|
| 1 aceptada | 21 | 7 | 28 |
| 3 una componente | 3 | 19 | 22 |
| total | 24 | 26 | 50 |

Acuerdo 40/50 (80 %), **κ = 0.60**; por categoría: `bd_ok` 20/27,
`bulbo_no_resuelto` 14/16, `disco_puro` 6/7. Ninguna de las seis tandas se aparta
de las demás (6/8 a 9/10, muestras demasiado chicas para distinguirlas). El
desacuerdo va en un sentido: en 7 de 10, aceptada a mano y una componente para
Claude.

**Revisiones dirigidas:**

| Lista | Revisadas | Resultado |
|---|---|---|
| Rechazadas (5) y dudosas (8) de Claude | 13 | las 5 rechazadas, confirmadas; de las 8 dudosas, 3 rechazadas, 3 aceptadas y 2 de una componente |
| Con `b` (20 al azar de 101) | 20 | clase cambiada en 5 y `b` quitada en 5: Claude marca barras de más en ~¼ de los casos |
| "Tipo elíptica" (B/T > 0.5, n_b ≥ 2, una componente para Claude; 15 al azar) | 15 | las 15, de una componente: coincide con Claude |

**Ajustes tras la revisión:**

- **Disco partido** (n_b < 1, R_e,b/R_e,d ≥ 0.5, B/T > 0.5): en todas las rondas
  la inspección humana los acepta y Claude los ponía de una componente. Se
  reclasificaron como aceptadas las 14 sin decisión humana, cada una con una
  nota. Copia previa en `inspeccion_visual_claude_antes_reclasificar.csv`.
- **`disco_puro` aceptadas:** a mano, "aceptada" en un `disco_puro` significa
  "el ajuste de una componente es bueno". En `visual_ok_final` se convierte en 3.
- **Grupo "tipo elíptica":** se probó un criterio físico (aceptar si el disco
  supera al bulbo por más de 0.6 mag a 2 R_e,d) para detectar S0 entre las de
  una componente de Claude. En las 15 nuevas acertó solo 9 (60 %), y no se
  aplicó. En las 30 galaxias del grupo revisadas a mano, Claude coincide en el
  90 %, así que el grupo queda sin cambios. Un corte en B/T (< 0.63) se descartó
  porque seleccionar por B/T, ligado a la masa del bulbo, puede sesgar la
  pendiente masa–tamaño. Valores del exceso del disco en `exceso_disco_claude.csv`.
- **Barras:** no se revisaron las 81 restantes. En el análisis se prueba la
  relación masa–tamaño de los bulbos con y sin las galaxias con `b`.

### Clasificación final (`--unir --final`)

```bash
python3 inspeccion_visual_bd.py --unir --final
```

Agrega a `galfitm_bd_filtrado.csv` y `componentes_masas.csv`:

| Columna | Contenido |
|---|---|
| `visual`, `visual_ok`, `visual_etiquetas` | solo las decisiones humanas (como `--unir`) |
| `visual_ok_final` | decisión humana; si no hay, la de Claude; si no, la automática (3) |
| `visual_fuente` | `humana`, `claude` o `auto` |
| `visual_etiquetas_final` | etiquetas con la misma prioridad (la barra aparece como `barra`) |

**El análisis usa `visual_ok_final`.** Con `--unir` sin `--final` se quitan las
tres columnas finales. Los `.bak` se escriben solo la primera vez (versiones
anteriores a cualquier unión).

`componentes_masas.csv` (1287 galaxias):

| `visual_fuente` | 0 rechazada | 1 aceptada | 3 una comp. | total |
|---|---|---|---|---|
| humana | 12 | 99 | 104 | 215 |
| claude | 0 | 459 | 508 | 967 |
| auto | 0 | 0 | 105 | 105 |
| total | 12 | 558 | 717 | 1287 |

`galfitm_bd_filtrado.csv` (4072 galaxias): 251 humanas, 1283 de Claude, 138
automáticas y 2400 sin clase (categorías que no se inspeccionan: `fallida`,
`invertida`, etc.); en total 735 aceptadas, 923 de una componente y 14 rechazadas.
Ninguna queda como dudosa.

Aceptadas con `b`: 74 de 558 en `componentes_masas.csv` (13 %; 51 de Claude, 23
humanas), que son las que se quitan en la prueba de robustez.

### Resumen para el paper

| Ronda | Galaxias | κ (4 clases) | Acuerdo |
|---|---|---|---|
| 1 (a ciegas, criterios generales) | 101 | 0.51 | 70 % |
| 2 (a ciegas, reglas finales) | 52 | 0.66 | 83 % |
| 3 (a ciegas, 50 al azar del paso 2) | 50 | 0.60 | 80 % |

De las 1287 galaxias de `componentes_masas.csv`, el 17 % tiene decisión humana,
el 75 % de Claude y el 8 % automática.

### Correcciones de código hechas durante la calibración

| Script | Problema | Arreglo |
|---|---|---|
| `inspeccion_visual_bd.py` | fallaba con galaxias de borde todo en ceros (isofotas) | respaldo de 3 líneas para estimar el ruido |
| `filtro_una_componente.py` | tomaba `RE_B_FWHM` como FWHM de la PSF: `r_disco_domina_arcsec` se buscaba desde ~0.14″ | usa `FWHM_PSF`. Ninguna categoría cambió (64 `sin_cruce` y 75 `disco_duplicado`, las mismas); cambió `r_disco_domina_arcsec` y, en 13 galaxias, si el cruce cae dentro de 1 FWHM |
| `inspeccion_visual_bd.py` | con dos ventanas abiertas, una podía sobrescribir lo que guardó la otra | cada decisión pasa por `guardar_decision`, que vuelve a leer el CSV y cambia solo la fila de esa galaxia; ya no se guarda la tabla entera al cerrar |
| `inspeccion_visual_bd.py` | `--unir` solo usaba las decisiones humanas | opción `--final` (ver *Clasificación final*); sin ella, la salida es idéntica a la anterior |
| `paso2_claude.py` | rearmar `inspeccion_visual_claude.csv` desde los lotes desharía la reclasificación del disco partido | se niega a rearmarlo |

---

## Aplicar a la muestra completa

Todo lo anterior se hizo con una **submuestra** (`galaxias_1_9` +
`galaxias_10_22`). Al sumar carpetas nuevas:

1. **Ajustes y cortes** (pasos 1 a 8): `python3 actualizar_todo.py --nproc 3`.
   Como algunos umbrales de `procesar_muestra.py` dependen de la población, unas
   pocas galaxias pueden entrar o salir de A + B. Las que salen conservan su
   clasificación visual en los CSV, pero `--unir` solo la pega a las galaxias que
   siguen en las tablas.
2. **Clasificación lógica** (pasos 9, 10 y 12, y el 11 si se usa):
   ```bash
   python3 leer_output_bd.py
   python3 filtro_bd.py
   python3 filtro_una_componente.py
   ```
   Los umbrales de estos dos filtros son fijos, así que las galaxias ya
   clasificadas no cambian de categoría al sumar otras.
3. **Ronda de control a ciegas**, antes de que Claude clasifique las nuevas:
   clasificar a mano 30–50 galaxias **de las carpetas nuevas** (`--orden azar`
   con otra semilla) y pedir a Claude la misma comparación. Las carpetas nuevas
   pueden venir de zonas de DP2 con otra profundidad o PSF, así que el κ de la
   submuestra no se puede dar por hecho.
4. **Clasificación de las nuevas por Claude**: solo las que no estén en
   `inspeccion_visual.csv` ni en `inspeccion_visual_claude.csv`, con las mismas
   reglas, agregando al mismo archivo.
5. **Control humano**: el mismo esquema (rechazadas, dudosas, `b` y 50 al azar)
   sobre las nuevas.
6. `python3 inspeccion_visual_bd.py --unir --final`. Repetirlo también si se
   vuelve a correr `filtro_bd.py` o `filtro_una_componente.py`.

Mensaje para Claude Code al sumar carpetas (la sesión nueva no recuerda las
anteriores, así que lleva todo el contexto):

```
You will act as my visual classifier for GalfitM bulge+disc decompositions, continuing
work done earlier on a subsample. Read README.md, section "Paso 14", first: it has the
classification codes, the general criteria, the final rules and the calibration
results. Apply exactly those rules.

New galaxies have been added (new galaxias_* folders). Steps:
1. CONTROL ROUND: the objectIds in analisis_galfitm/bd/filtro/inspeccion_visual.csv
   from the new folders are my blind control set. Classify them from --limpia figures
   WITHOUT reading my visual, visual_ok or visual_etiquetas columns, then report the
   confusion matrix, Cohen's kappa (all classes and accepted vs. rest) and the
   disagreements with your reasoning. STOP and show me this.
2. If I approve, classify only the galaxies that are in neither inspeccion_visual.csv
   nor inspeccion_visual_claude.csv (default categories, without --incluir-auto), in
   saved batches of ~50, appending to inspeccion_visual_claude.csv with
   clasificador=claude and a short note. Do NOT modify inspeccion_visual.csv.
3. At the end, give me as separate CSVs: your rejected and uncertain galaxies, all
   galaxies with flag b, and a random sample of 50 of the rest (fixed seed), plus a
   check of 20 random auto-classified galaxies from the new folders.
```

Para el paper: el prompt completo en inglés va en un apéndice, y en métodos se
describe el procedimiento (criterios, calibración a ciegas con el κ obtenido y
fracción revisada a mano).

---

## Uso diario: `actualizar_todo.py`

Desde `analisis_conjunto/`. Recorre las carpetas `galaxias_*` (junto a esta
carpeta o un nivel arriba) y corre solo lo pendiente:

1. single-Sérsic de galaxias sin salida (no reintenta las fallidas) y
   `leer_output_galfitm.py` si hay salidas nuevas;
2. `combinar_carpetas.py` y `procesar_muestra.py` (cachés: solo mide lo nuevo);
3. statmorph y bulbo + disco **solo** para las A + B que no los tienen;
4. une los resultados de statmorph.

La clasificación bulbo + disco (pasos 9 a 14) **no** la corre: después de
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

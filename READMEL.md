# procesar_muestra.py y corregir_ext_kcor.py

Selección de una muestra confiable de galaxias a partir de ajustes multibanda de
**GalfitM** sobre imágenes de **LSST**, validados contra el catálogo de LSST
(Sérsic de MultiProFit) y contra las propias imágenes.

El script reúne en un solo archivo tres etapas del análisis:

| Paso | Qué hace | Entrada | Salida principal |
|------|----------|---------|------------------|
| **2** | Cortes de calidad de GalfitM | `galfitm_resultados.csv` | `analisis_galfitm/galfitm_limpio.csv` |
| **3** | Cruce (*crossmatch*) con el catálogo LSST | `galfitm_limpio.csv` + parquet LSST | `analisis_galfitm/crossmatch/galfitm_x_lsst.csv` |
| **4** | Clasificación de confiabilidad A/B/C/D | `galfitm_x_lsst.csv` + FITS de GalfitM | `analisis_galfitm/muestra/muestra_clasificada.csv` |

Después, `corregir_ext_kcor.py` corrige las magnitudes de esa muestra por
extinción galáctica y corrección K, y calcula magnitudes absolutas, colores en
reposo y masa estelar (ver la [sección correspondiente](#corregir_ext_kcorpy--extinción-corrección-k-y-masa-estelar)).

---

## Contexto: el pipeline completo

```
correr_galfitm.sh          ->  ajusta GalfitM a cada galaxia (*_galfitm_out.fits)
leer_output_galfitm.py     ->  paso 1: lee los FITS y crea galfitm_resultados.csv
procesar_muestra.py        ->  pasos 2, 3 y 4  (este script)
corregir_ext_kcor.py       ->  extinción galáctica, corrección K, magnitudes absolutas, masa
clasificar_galaxias.py     ->  clasificación color + índice de Sérsic
```

---

## Requisitos

- Python ≥ 3.9
- `numpy`, `pandas`, `matplotlib`, `astropy`, `pyarrow` (para leer el parquet)

```bash
pip install numpy pandas matplotlib astropy pyarrow
```

Para `corregir_ext_kcor.py` además:

- `calc_kcor_fun.py`: módulo de corrección K de Chilingarian et al.
  ([kcor.sai.msu.ru](http://kcor.sai.msu.ru)), en la misma carpeta que el script.
- `dustmaps` (opcional): solo si se recalcula E(B−V) con el mapa SFD
  (`EBV_FUENTE = 'sfd'`).

## Estructura de carpetas esperada

```
mi_muestra/
├── procesar_muestra.py
├── corregir_ext_kcor.py
├── calc_kcor_fun.py                # módulo de Chilingarian (para la corrección K)
├── galfitm_resultados.csv          # salida de leer_output_galfitm.py
├── z01_lsst_1.parquet              # catálogo LSST (ruta configurable)
└── galfitm_inputs/
    └── <objectId>/
        ├── <objectId>_galfitm_out.fits   # INPUT_b, MODEL_b, RESIDUAL_b por banda
        ├── <objectId>_mask.fits
        └── <objectId>_sig_<banda>.fits
```

Todas las rutas son relativas a la ubicación del script, así que basta con
copiarlo en la carpeta de cada muestra.

## Uso

```bash
python3 procesar_muestra.py          # corre los pasos 2, 3 y 4
python3 procesar_muestra.py 4        # solo el paso 4 (p. ej., tras cambiar un umbral)
python3 procesar_muestra.py 3 4      # pasos 3 y 4
```

Antes de empezar, el script comprueba que existan el parquet y la salida de cada
paso anterior; si falta algo, indica qué archivo y se detiene.

## Configuración general

Al inicio del script:

| Variable | Por defecto | Descripción |
|----------|-------------|-------------|
| `LSST_PARQUET` | `z01_lsst_1.parquet` | Ruta al catálogo LSST |
| `REVISAR_IMAGENES_FITS` | `True` | Paso 2: abre los FITS para detectar bandas vacías o con fondo malo |
| `RECALCULAR_METRICAS` | `False` | Paso 4: fuerza a rehacer las métricas de imagen |

Los umbrales de cada corte están al inicio de cada paso (buscar `CONFIG` en el código).

---

## Paso 2 — Cortes de calidad de GalfitM

Aplica cortes en dos niveles.

**Por galaxia** (se descarta la galaxia completa):

| Corte | Valor | Motivo |
|-------|-------|--------|
| χ²/ν mínimo | `CHI2_MIN = 0.5` | χ²/ν ≈ 0 indica un problema (sigma mal estimada, banda vacía) |
| χ²/ν máximo | `CHI2_MAX = 5.0` | el modelo no reproduce la galaxia |

**Por banda** (solo esa banda pasa a NaN; el resto de la galaxia se conserva):

| Motivo | Condición |
|--------|-----------|
| `flag_galfitm` | GalfitM marcó algún parámetro con `*` |
| `sin_error` | error de magnitud vacío o ≤ 0 |
| `mag_fuera` | magnitud fuera de [10, 25] |
| `re_fuera` | R_e fuera de [0.5, 150] px |
| `n_limite` | n fuera de [0.21, 7.9] (pegado a los límites de GalfitM) |
| `ar_bajo` | b/a < 0.1 |
| `color_outlier` | la magnitud de la banda se aleja > 5σ (y > 1 mag) del color típico de la población |
| `imagen_vacia` | > 30 % de píxeles = 0 o NaN en la imagen INPUT |
| `fondo_malo` | variación del fondo a gran escala > 3 veces el ruido (gradientes, imágenes sin ruido) |

Las dos últimas solo se aplican con `REVISAR_IMAGENES_FITS = True`. El ruido se
mide con diferencias entre píxeles vecinos (insensible a gradientes) y el fondo
con medianas en bloques de 20×20 píxeles.

**Salidas:** `galfitm_limpio.csv`, `galaxias_a_revisar.csv`, `calidad_imagenes.csv`
y las figuras 01–07.

---

## Paso 3 — Crossmatch con LSST

- Cruza por `objectId` (como entero de 64 bits, para no perder precisión).
- Si un objeto LSST tiene varios espectros, conserva uno: prioriza `f_z = 1` y,
  entre esos, la menor distancia `_dist_arcsec`.
- Elimina los valores centinela de `cModelMag` (≈ 80).
- Convierte `{b}_sersicFlux` (nJy) a magnitud AB: `m = −2.5 log10(F) + 31.4`.
  Es el mismo zeropoint de los ajustes de GalfitM, así que ambas magnitudes son
  directamente comparables.
- Compara:
  - magnitud de GalfitM con `cModelMag` y con la magnitud Sérsic de LSST, por banda;
  - R_e, n y b/a de GalfitM (banda i) con MultiProFit (solo objetos sin
    `shape_flag` ni flags de fallo).

**Salidas:** `galfitm_x_lsst.csv`, `resumen_magnitudes.csv`, `sin_cruce.csv`
(si hay galaxias sin pareja) y las figuras 08–11.

---

## Paso 4 — Muestra confiable (categorías A/B/C/D)

**Los parámetros de LSST no se toman como la verdad.** MultiProFit es otro ajuste
paramétrico (un solo n para todas las bandas, n limitado a [0.5, 6], afectado por
el *deblender*), así que también puede equivocarse. LSST se usa como una
**segunda medida independiente**.

### Lógica

1. **Si GalfitM y LSST coinciden**, dos códigos distintos llegan al mismo resultado.
2. **Si no coinciden**, al menos uno se equivoca. Para decidir, se usa un árbitro
   independiente de ambos modelos: **la imagen**. En una elipse de 2 R_e se mide:
   - el **flujo de apertura** sumando directamente los píxeles (sin modelo);
   - el cociente **flujo del modelo de GalfitM / flujo de la imagen**;
   - el **RFF** (*residual flux fraction*, Hoyos et al. 2011): flujo del residuo
     por encima del esperado por el ruido.

### Acuerdo con LSST

Las diferencias GalfitM − LSST se miden **respecto al offset mediano** de la
muestra, para que un sesgo sistemático no elimine galaxias:

| Parámetro | Tolerancia |
|-----------|------------|
| Magnitud | ± 0.5 mag |
| log R_e | ± 0.15 dex |
| log n | ± 0.15 dex (no se compara si el n de LSST está en sus límites, 0.5 o 6) |
| b/a | ± 0.10 |

La comparación usa la banda i; si falta, r, z, g o y, en ese orden.

### Requisitos previos (si no se cumplen → D)

- Pasa los cortes de galaxia del paso 2 y tiene al menos 3 bandas buenas.
- Tiene buenas las **bandas obligatorias** `['g', 'r', 'i']`:
  g y r para la clasificación (g − r, n_r), i para la masa estelar de
  Taylor et al. (2011), que usa (g − i) y M_i.
- Espectro confiable: `f_z = 1`, clase `GALAXY`, distancia del cruce ≤ 1″.

### Categorías

| Categoría | Condición | Interpretación |
|-----------|-----------|----------------|
| **A** | GalfitM y LSST coinciden | Alta confianza: dos códigos independientes concuerdan |
| **B** | No coinciden (o LSST no es comparable), pero GalfitM reproduce la imagen: RFF < 0.10 y flujo modelo/imagen entre 0.9 y 1.1 | Confiable; la diferencia se atribuye a LSST o al método |
| **C** | No coinciden y GalfitM tampoco reproduce bien la imagen | Revisar a ojo |
| **D** | No cumple los requisitos previos | Descartada |

La columna `motivo_categoria` explica por qué cada galaxia está en su categoría,
y `problemas_lsst` anota los problemas internos del ajuste de LSST
(`shape_flag`, galaxia fragmentada por el *deblender*, n en el límite, χ² alto).
Estos últimos son informativos y no deciden la categoría.

### Métricas de imagen en caché

La primera ejecución abre todos los FITS (tarda unos minutos) y guarda los
resultados en `metricas_imagen.csv`; las siguientes los reutilizan. Si las
métricas guardadas no cubren las galaxias actuales (por ejemplo, porque son de
otra muestra), el script lo detecta y las recalcula automáticamente.

**Salidas:** `muestra_clasificada.csv` (todas las galaxias), `muestra_A.csv`,
`muestra_AB.csv`, `metricas_imagen.csv` y las figuras 12–14.

---

## corregir_ext_kcor.py — Extinción, corrección K y masa estelar

Se ejecuta después de `procesar_muestra.py`:

```bash
python3 corregir_ext_kcor.py
```

Usa las categorías **A + B** de `muestra_clasificada.csv` (si no existe, usa
`galfitm_x_lsst.csv` con `f_z = 1`). Aplica el mismo método que el pipeline de
S-PLUS, para que los resultados de ambas muestras sean comparables.

### 1. Extinción galáctica

```
A_V = R_V · E(B−V)        R_V = 3.1
A_λ = CCM89(λ_eff, A_V)   Cardelli, Clayton & Mathis (1989)
m_ext = m_GalfitM − A_λ
```

- **E(B−V):** por defecto la columna `ebv` del catálogo LSST
  (`EBV_FUENTE = 'columna'`); con `EBV_FUENTE = 'sfd'` se recalcula con el mapa
  SFD mediante `dustmaps`. Conviene verificar que la columna de LSST provenga de SFD.
- **Escala de SFD:** `SFD_ESCALA = 1.0` usa SFD original (igual que S-PLUS);
  `0.86` aplica la recalibración de Schlafly & Finkbeiner (2011). Debe usarse el
  mismo valor en todas las muestras.
- **CCM89 completo:** incluye la rama infrarroja (x < 1.1 μm⁻¹), necesaria para
  la banda y. En u–z coincide con la función óptica del pipeline de S-PLUS.

Longitudes de onda efectivas LSST (Ivezić et al. 2019) y coeficientes resultantes:

| Banda | λ_eff [Å] | A_λ / E(B−V) |
|-------|-----------|--------------|
| u | 3671 | 4.81 |
| g | 4827 | 3.64 |
| r | 6223 | 2.70 |
| i | 7546 | 2.06 |
| z | 8691 | 1.58 |
| y | 9712 | 1.31 |

### 2. Corrección K

Con el método de Chilingarian, Melchior & Zolotukhin (2010), usando colores ya
corregidos por extinción. Para cada banda se prueba la lista de colores en orden
y se usa el primero disponible en cada galaxia:

| Banda | Colores (en orden de preferencia) |
|-------|-----------------------------------|
| u | u − r, u − i, u − z |
| g | g − r, g − i, g − z |
| r | g − r, u − r |
| i | g − i, u − i |
| z | g − z, r − z, u − z |

La columna `kcor_color_<banda>` registra qué color se usó en cada galaxia.
La banda **y** no tiene polinomio en `calc_kcor`: queda corregida solo por
extinción (sin M_y).

### 3. Magnitudes absolutas y colores en reposo

```
M = m_ext − K(z, color) − DM(z)        cosmología plana, H0 = 70, Ωm = 0.3
(g − r)_0 = M_g − M_r      (u − r)_0 = M_u − M_r      (g − i)_0 = M_g − M_i
```

Se usa el redshift espectroscópico `z`; para z < 0.001 la distancia queda en NaN.

### 4. Masa estelar

Taylor et al. (2011), con color y magnitud **en reposo**:

```
log10(M*/M☉) = 1.15 + 0.70 (g − i)_0 − 0.40 M_i
```

Como la banda i es obligatoria en `procesar_muestra.py`, todas las galaxias A + B
tienen masa.

### Salidas

En `analisis_galfitm/fotometria/`:

| Archivo | Contenido |
|---------|-----------|
| `galaxias_ext_kcor.csv` | Por banda: `MAG_` (GalfitM), `A_` (extinción), `mag_ext_`, `kcor_`, `kcor_color_`, `M_`. Además `EBV_usado`, `DM`, colores observados (`*_obs`) y en reposo (`*_0`), y `logM_taylor` |
| `17_kcorrecciones.png` | Corrección K vs z para cada banda |
| `18_color_observado_vs_reposo.png` | Color vs z antes y después de la corrección K; la línea roja sigue a las galaxias con n_r ≥ 2.5. Si la corrección funciona, la tendencia con z desaparece o se reduce |

### Advertencia

Los polinomios de `calc_kcor` están calibrados para los filtros de **SDSS**. Los
de LSST son similares y a z < 0.1 la diferencia es pequeña, pero se trata de una
aproximación que conviene mencionar.

---

## Figuras

| Figura | Contenido |
|--------|-----------|
| `01_chi2.png` | Distribución de χ²/ν con los cortes |
| `02_histogramas_por_banda.png` | Magnitud, R_e, n y b/a por banda, antes y después de los cortes |
| `03_dependencia_lambda.png` | Variación de R_e, n y b/a con la longitud de onda |
| `04_colores.png` | Color-magnitud, color-color y color vs n |
| `05_n_re_ar.png` | n vs R_e y b/a vs n |
| `06_errores.png` | Errores formales de GalfitM vs magnitud |
| `07_calidad_imagenes.png` | Métricas de fondo y fracción de píxeles vacíos (para calibrar los cortes) |
| `08_mag_vs_cmodel.png` | GalfitM vs `cModelMag` de LSST, por banda |
| `09_mag_vs_sersic_lsst.png` | GalfitM vs Sérsic de LSST, por banda |
| `10_estructura_vs_lsst.png` | R_e, n y b/a: GalfitM vs MultiProFit |
| `11_dmag_vs_n.png` | Diferencia de magnitud vs índice de Sérsic |
| `12_consistencia_lsst.png` | Diferencias GalfitM − LSST con las tolerancias |
| `13_arbitro_imagen.png` | RFF, cociente de flujo y comparación con la fotometría de apertura |
| `14_categorias.png` | Propiedades de las galaxias en cada categoría |

Las figuras 17 y 18 las genera `corregir_ext_kcor.py` (ver su sección).
Los títulos y ejes de todas las figuras están en inglés.

---

## Advertencias

- **Sesgo por bandas obligatorias.** Exigir la banda i elimina ~10 % de la
  muestra. Esas galaxias son algo más rojas y esferoidales que el resto (les
  falta i por la cobertura de LSST, no por mala calidad). Conviene mencionarlo
  en la descripción de la muestra.
- **Sesgo de la categoría A.** Las galaxias grandes, brillantes y cercanas suelen
  quedar fuera de A porque el *deblender* de LSST las fragmenta. Para estudios
  de estructura se recomienda usar A + B.
- **Factor de escala en R_e.** El R_e de GalfitM es sistemáticamente ~1.46 veces
  el de MultiProFit, sin dependencia con b/a ni con n. Es una diferencia de
  escala o de definición que debe resolverse antes de publicar.
- **Errores formales.** Los errores de GalfitM están subestimados; el header
  guarda solo 4 decimales (piso de 10⁻⁴ mag).
- **Columnas de LSST.** Las columnas informativas (`shape_flag`,
  `deblend_blendNChild`, `sersic_chi2_reduced`) son opcionales; si faltan en el
  parquet se usan valores neutros. Las columnas del espectro (`f_z`, `class`,
  `_dist_arcsec`) y del Sérsic de LSST sí son necesarias.

## Referencias

- Cardelli, J. A., Clayton, G. C. & Mathis, J. S. (1989), ApJ, 345, 245 — ley de extinción CCM89
- Chilingarian, I. V., Melchior, A.-L. & Zolotukhin, I. Yu. (2010), MNRAS, 405, 1409 — corrección K
- Chilingarian, I. V. & Zolotukhin, I. Yu. (2012), MNRAS, 419, 1727 — extensión de la corrección K
- Häußler, B. et al. (2013), MNRAS, 430, 330 — GalfitM / MegaMorph
- Hoyos, C. et al. (2011), MNRAS, 411, 2439 — *Residual flux fraction*
- Ivezić, Ž. et al. (2019), ApJ, 873, 111 — LSST (filtros)
- Schlafly, E. F. & Finkbeiner, D. P. (2011), ApJ, 737, 103 — recalibración de SFD
- Schlegel, D. J., Finkbeiner, D. P. & Davis, M. (1998), ApJ, 500, 525 — mapa de polvo SFD
- Taylor, E. N. et al. (2011), MNRAS, 418, 1587 — masa estelar a partir de (g − i)

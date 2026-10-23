# deku2d — esquemas 2D de interacciones ligando–receptor

Genera figuras 2D de calidad de artículo de la unión de un ligando a su receptor (pose de docking o fotograma de
dinámica molecular). Las interacciones las calcula **PLIP**; el ligando se dibuja con **RDKit** y todo el diseño y el
dibujo se hacen con **matplotlib**. El estilo reproduce `estilo.png`:

| Elemento | Representación |
|---|---|
| Ligando | trazo negro fino; solo los heteroátomos llevan letra |
| Contacto hidrofóbico / van der Waals | halo amarillo sobre el átomo o anillo, 2 arcos naranjas hacia el residuo y línea punteada gris clara caja→arcos |
| Residuo | caja gris con nombre en mayúsculas y cadena (`ASP12A`) |
| Enlace de H, dador en la proteína | línea punteada **roja** con flecha triangular abierta hacia el ligando |
| Enlace de H, dador en el ligando | línea punteada **verde** con flecha hacia el residuo |
| Puente salino | línea punteada morada |
| π–π / π–catión | línea punteada cian / naranja con punto en el centroide |
| Distancia (Å) | solo en enlaces de H y puentes salinos, con el color de la línea (`--no-dist` la quita) |

## Requisitos

* `python3` con **rdkit**, **matplotlib** y **numpy** (en este equipo ya están en el `python3` del sistema).
* **PLIP 2.x** solo si se parte de un PDB: se llama como programa externo, por defecto
  `/home/alecastro/anaconda3/envs/plip_env/bin/plip` (otra ruta: `--plip-exe` o la variable `DEKU2D_PLIP`).
  Con un `report.xml` ya calculado no hace falta PLIP.

## Uso

Sin instalar, desde cualquier carpeta:

```bash
/home/alecastro/Documentos/proyecto_2Deku_scheme/deku2d.sh --help
# opcional: tenerlo siempre disponible como «deku2d»
ln -s /home/alecastro/Documentos/proyecto_2Deku_scheme/bin/deku2d ~/.local/bin/deku2d
# o instalar:  pip install -e /home/alecastro/Documentos/proyecto_2Deku_scheme   (crea la orden «deku2d»)
# o desde la carpeta del proyecto:  python3 -m deku2d --help
```

### Un esquema (`draw`)

```bash
# A partir del complejo PDB (ejecuta PLIP; el SMILES define enlaces y cargas)
deku2d draw --pdb complejo.pdb --ligand-resname UNK --smiles "OCc1[nH]c(O)cc1..." -o figura.png

# A partir de un report.xml de PLIP ya calculado (usa el plipfixed*.pdb que está junto al XML)
deku2d draw --plip-xml ejemplos/mini/mol3_report.xml --ligand-resname UNK --smiles "..." -o figura.png

# Ligando como SDF (se asocia a los átomos del PDB por coordenadas) y varios formatos a la vez
deku2d draw --pdb complejo.pdb --ligand-resname 6IC --sdf ligando.sdf --format png,svg,pdf -o figura.png
```

### Panel con varios sistemas (`panel`) y lote (`batch`)

```bash
deku2d panel --config ejemplos/panel_docking_cuatro.json -o fig_panel.png          # a), b), c), d) + una leyenda
deku2d batch --config ejemplos/batch_ejemplo.json --outdir salida                  # una figura por sistema (+ panel)
```

JSON de panel/lote (las rutas relativas se resuelven respecto al JSON):

```json
{
  "ncols": 2, "label_style": "a)", "formats": ["png", "svg"],
  "style": {"lang": "en", "dist": true},
  "systems": [
    {"title": "mol3 (docking)", "name": "mol3", "plip_xml": "mol3_report.xml", "ligand_resname": "UNK", "smiles": "..."},
    {"title": "mol3 (MD 200 ns)", "name": "mol3_md", "pdb": "mol3_rep.pdb", "ligand_resname": "UNK", "smiles": "..."}
  ],
  "panel": {"output": "panel", "ncols": 2}
}
```

Claves de cada sistema: `ligand_resname` (obligatoria), `pdb`, `smiles`, `sdf`, `plip_xml`, `chain`, `resnum`, `charge`, `title`, `name`, `workdir`, `plip_exe`.

### Opciones principales

`--lang en|es` (leyenda), `--dist/--no-dist`, `--fontsize` (17), `--ppu` (34; menor = letra relativamente mayor),
`--dpi` (300), `--format png,svg,pdf`, `--min-width`, `--seed`, `--restarts`/`--iterations` (búsqueda de colocación),
`--no-legend`, `--title`, `--style-json` (cualquier clave de `deku2d.config.Style`), `--legacy` (ver abajo), `-q`.
El `-o` fija el nombre base; con `--format png,svg` se escriben `figura.png` y `figura.svg`.
Prioridad: valores por defecto < `style` del JSON < `--style-json` < opciones de la línea de órdenes.

## Datos que necesita cada sistema

1. **Estructura del complejo** (PDB) donde el ligando es un residuo `HETATM` (o un `report.xml` de PLIP + su `plipfixed*.pdb`).
   Los índices de átomo del XML son números de serie del PDB **corregido por PLIP**, no los del PDB original: por eso, si se
   da `--plip-xml`, se usa el `plipfixed*.pdb` de su carpeta. Si no existe, se usa `--pdb` y sus series deben coincidir.
2. **Código de residuo del ligando** (`--ligand-resname`, p. ej. `UNK`, `6IC`). Si hay varias copias, `--chain`/`--resnum`.
3. **Estructura química del ligando**, para dibujar enlaces dobles, aromaticidad y cargas:
   * **SMILES plantilla** (`--smiles`, recomendado): mismos átomos pesados que el PDB (los hidrógenos del SMILES son
     implícitos); el orden no importa. Los enlaces y las **cargas formales salen del SMILES**: escríbalo con el estado de
     protonación y el tautómero usados en el cálculo (p. ej. carboxilato `C(=O)[O-]`, amonio `[NH+]`). La estereoquímica se
     ignora (el esquema es plano).
   * **SDF/MOL** (`--sdf`): con enlaces y cargas; se asocia por elemento y coordenadas (o por orden si no coinciden).
   * **Nada**: se intentan deducir los enlaces de la geometría 3D (necesita hidrógenos) y, si falla, se dibujan enlaces
     simples con un aviso. Úselo solo para una vista previa.
4. Si el SMILES no encaja con el PDB (distinto número de átomos pesados o conectividad), la app se detiene con un mensaje
   que indica la diferencia.

### Docking y dinámica molecular

* **Docking**: la pose acoplada como PDB del complejo (proteína + ligando).
* **Dinámica**: exporte un fotograma representativo (p. ej. el medoide del grupo más poblado de la parte estable de la
  trayectoria) como PDB con proteína + ligando (el agua y los iones pueden omitirse). Dibuje con el mismo SMILES y
  `--ligand-resname`. Las interacciones corresponden a **ese** fotograma: es una instantánea, no una ocupación.

## Qué hace el diseño

Las cajas de residuos se colocan sobre una envolvente que sigue el contorno de la molécula (dos capas) mediante una
búsqueda local con reinicios que minimiza: longitud de las líneas, cruces con enlaces del ligando y entre líneas,
líneas que atraviesan cajas ajenas, cajas demasiado próximas y distancia a los arcos hidrofóbicos; las cajas nunca se
solapan. **Mejoras respecto al script original** (activas por defecto, `--legacy` las desactiva):

* las etiquetas de distancia se colocan donde menos estorban (a un lado u otro de la línea, evitando otras etiquetas, cajas,
  líneas y enlaces): ya no se pisan cuando varios enlaces llegan al mismo átomo;
* las líneas grises caja→arcos y las líneas polares se penalizan si cruzan el ligando o las letras de sus heteroátomos;
* el coste se evalúa con tablas y operaciones vectorizadas de numpy: ≈ 4 s por sistema (30 reinicios × 2500 movimientos)
  frente a ≈ 1 min del script original; con `--legacy` el resultado es **idéntico píxel a píxel** al del script original.

## Pruebas

```bash
cd /home/alecastro/Documentos/proyecto_2Deku_scheme
python3 -m unittest discover -s tests          # también funciona con pytest, si está instalado
```

Cubren el análisis del XML, el mapeo de átomos, `draw` (PNG/SVG/PDF), los paneles, la CLI y la ejecución real de PLIP
(se omite si PLIP no está). `ejemplos/mini/` contiene el complejo mol3, su `report.xml` y un panel de ejemplo;
`ejemplos/datos/` tiene sistemas adicionales.

## Limitaciones

* Solo se dibujan enlaces de H, contactos hidrofóbicos, puentes salinos, π–π y π–catión. Puentes de agua, enlaces de halógeno
  e interacciones con metales (PLIP) se ignoran con un aviso.
* El diseño es plano y automático: puede haber líneas que cruzan el ligando cuando varios residuos interactúan con átomos
  opuestos, y las etiquetas de dos enlaces muy cortos al mismo átomo pueden quedar juntas. Hay semilla (`--seed`) y
  parámetros de búsqueda para probar otra colocación.
* Las interacciones son las que detecta PLIP con sus umbrales por defecto sobre una estructura estática; no sustituyen un
  análisis de ocupación en la dinámica ni una estimación de energía.
* Los números de residuo y la cadena son los del PDB de entrada.
* Si dos ligandos tienen el mismo código, indique `--chain` o `--resnum`.

## Estructura

```
deku2d/            config (estilo, textos en/es) · plip_io · ligand · layout · render · panel · pipeline · cli
tests/             unittest (compatibles con pytest)
ejemplos/mini/     ejemplo mínimo (mol3) y panel.json
ejemplos/*.json    paneles/lotes de ejemplo (usan ejemplos/datos)
deku2d.sh, bin/deku2d   lanzador
```

Al usar PLIP en un artículo, cítelo (Adasme et al., *Nucleic Acids Res.* 2021, 49, W530–W534) y mencione RDKit.

## Figuras del manuscrito (KRAS G12D, Switch-II)

Los paneles de las figuras 12, 22 y 23 del manuscrito se generan así (los datos están en `ejemplos/datos/`):

```bash
P=/home/alecastro/Documentos/proyecto_2Deku_scheme
$P/deku2d.sh panel --config $P/ejemplos/panel_docking_cuatro.json        -o fig34_plip2d_docking.png      # poses de docking (4)
$P/deku2d.sh panel --config $P/ejemplos/panel_md_representativo_cuatro.json -o fig35_plip2d_md_rep.png    # fotograma representativo (4 candidatos)
$P/deku2d.sh draw --plip-xml $P/ejemplos/datos/md_rep/MRTX1133/report.xml \
  --pdb $P/ejemplos/datos/md_rep/MRTX1133/frame_rep_protonated.pdb --ligand-resname 6IC \
  --smiles "$(python3 -c "import json;print(json.load(open('$P/ejemplos/datos/ligandos.json'))['MRTX1133']['smiles'])")" \
  --as-salt ASP12 -o fig36_plip2d_md_mrtx1133.png                                                                     # referencia MRTX1133
```

Nota sobre números de serie: si no hay `plipfixed*.pdb` junto al `report.xml` y los números de serie del informe no coinciden con el
PDB (PLIP puede renumerar), la aplicación deduce el desplazamiento comparando las coordenadas `<ligcoo>` del informe con las del
PDB (en MRTX1133 era −1) y avisa. Las interacciones de los tipos puente salino, π–π y π–catión usan solo los átomos del ligando
(`<lig_idx_list>`).

## Opciones añadidas (2026-10-20)

- Letra de los átomos del ligando: por defecto `--fontsize + 2` (`--atom-fontsize` la fija a mano).
- Relleno de las cajas por tipo de residuo (`--tint`), borde grueso en residuos ausentes de una estructura de referencia
  (`compare_xml` y `compare_label` por sistema en el JSON; `--compare-xml` en `draw`), paneles con celdas uniformes (`--uniform`),
  distancias hidrofóbicas (`--vdw-dist`) y cajas repetidas numeradas (THR58A·1, THR58A·2).
- `--as-salt ASP12` (o `"as_salt": ["ASP12"]` en el JSON): dibuja como puente salino los enlaces de H de ese residuo. Es una
  reclasificación manual (PLIP los informa como enlaces de H) y la aplicación lo avisa.

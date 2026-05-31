# FITOAPP — Pipeline dei dati reali

Come l'app passa da **una scansione UAV multispettrale pubblica** ai **focolai di
irrorazione**, al **tragitto di guida** e alle **statistiche** mostrate a schermo.

Tutto il processing è in un solo script riproducibile: [`tools/process_dataset.py`](tools/process_dataset.py).
Gli asset finali usati dall'app stanno in [`data/`](data/).

```
 Zenodo (dataset reale)                  tools/process_dataset.py                 data/  ──►  index.html
 ─────────────────────                   ────────────────────────                 ─────       ─────────
 datasheet.csv  (viti + Esca)  ──┐   1. lettura + filtro vigneto B7        ┌─►  targets.geojson   focolai sulla mappa
 ortho_b7.tif   (multispettrale)─┤   2. focolai (clustering viti malate)   ├─►  route.geojson     percorso del trattore
 trunk_*.zip    (.prj = CRS)   ──┘   3. tragitto (serpentina sui filari)   ├─►  field.png         overlay drone RGB
                                     4. NDVI + overlay (RGB / NDVI)         ├─►  ndvi.png          overlay NDVI
                                     5. valutazioni / statistiche           ├─►  meta.json         centro, bounds, stat
                                                                            └─►  data.js           bundle (anche file://)
```

---

## 0. Da dove vengono i dati

**Dataset**: *Precision viticulture dataset — Northern Spain, July 2022* —
Vélez, Ariza-Sentís & Valente (2024), **Zenodo [10.5281/zenodo.10362568](https://doi.org/10.5281/zenodo.10362568)**, licenza **CC-BY-4.0**.

Rilievo con drone **DJI Matrice 210 RTK + sensore multispettrale Micasense Altum**
(5 bande: Blu, Verde, Rosso, Red-Edge, NIR) su un vigneto in Galizia (Spagna).
Il dataset contiene due vigneti, **B7** e **B9**. Si usa **solo B7** perché
l'ortomosaico scaricato copre quel blocco; le viti di B9 cadrebbero fuori dall'overlay.

File scaricati (in `raw/`, **non** versionati perché pesanti):
| File | Cos'è | Uso |
|---|---|---|
| `datasheet.csv` | una riga per vite: coordinate UTM + colonna `Esca` (YES/NO) | **focolai** |
| `..._B7_..._MS_PRO.tif` (~3.2 GB) | ortomosaico multispettrale georeferenziato | **NDVI + overlay** |
| `trunk_locations_VineyardB7.zip` | shapefile (serve solo il `.prj` per leggere il CRS) | sistema di riferimento |

> Il "risultato del processing" che l'app sfrutta è già nel `datasheet.csv`: ogni
> pianta è **georeferenziata e diagnosticata** (malattia Esca sì/no). Da lì si parte.

---

## 1. Lettura delle viti e filtro

`load_vines()`:
1. Legge `datasheet.csv`, tenendo **solo le righe del vigneto B7** (`Vineyard == "B7"`).
2. Per ogni vite prende le coordinate UTM `X,Y` e lo stato `Esca`.
   Una vite è **malata** se `Esca` inizia con `YES`.

Risultato: **133 viti B7, di cui 38 malate**.

Sistemi di riferimento coinvolti (gestiti con `pyproj`/`rasterio`):
- viti → **EPSG:25829** (ETRS89 / UTM 29N, dal `.prj`)
- ortomosaico → **EPSG:32629** (WGS84 / UTM 29N)
- mappa Leaflet → **EPSG:4326** (lat/lng)

ETRS89 e WGS84 in Europa differiscono di <1 m: l'allineamento resta sub-metrico.

---

## 2. Come si trovano i focolai

Le **viti malate** sono i punti da irrorare, ma mostrarne 38 separate sarebbe
illeggibile: si raggruppano in **focolai** (cluster di infezione vicini).

`cluster()` + `foci_from_clusters()`:
1. **Clustering** con *union-find* (single linkage): due viti malate finiscono nello
   stesso focolaio se distano ≤ **`LINK_DIST = 6.5 m`**. Catene di viti contigue
   diventano un unico focolaio.
2. Per ogni focolaio:
   - **centro** = media delle coordinate delle viti del gruppo;
   - **raggio** = distanza massima centro→vite **+ `FOC_MARGIN = 2.0 m`**, limitato a
     **`RAD_CAP = 14 m`** (evita focolai allungati/giganti dovuti all'effetto catena);
   - **n. viti** del gruppo.
3. Il centro viene **riproiettato in lat/lng** (WGS84) per Leaflet.

Risultato: **11 focolai**. I valori di `LINK_DIST`/`RAD_CAP` sono stati scelti
cercando un compromesso tra leggibilità (pochi marker) e fedeltà (focolai compatti).

→ output **`data/targets.geojson`** (centro, raggio, n. viti). Nell'app questi sono i
`FOCI`: aloni arancioni sulla mappa e zone in cui la guida dà **SPRUZZA**.

---

## 3. Come si traccia il tragitto

Un trattore non si muove a caso: percorre i **filari** avanti e indietro
(*boustrophedon* / serpentina). Il problema è capire **come sono orientati i filari**.

`serpentine_route()`:
1. **Direzione dei filari via PCA** — si calcola la matrice di covarianza delle
   posizioni delle viti e se ne ricava l'asse principale: è la direzione lungo cui le
   viti sono più "allungate", cioè il **filare**. L'asse perpendicolare è la direzione
   *tra* i filari.
   ```
   θ = ½ · atan2(2·Σxy , Σxx − Σyy)     # angolo dell'asse principale
   ```
2. Si proiettano tutte le viti su questi due assi (u = lungo il filare, v = tra filari)
   per ottenere l'estensione reale del campo.
3. Si generano **5 passate**: per ciascuna due estremi (inizio/fine filare), con
   **direzione alternata** (serpentina). Un piccolo `pad` allarga le passate oltre le
   viti di bordo.
4. Gli estremi vengono riconvertiti in lat/lng.

→ output **`data/route.geojson`** (i *waypoint*). Nell'app la funzione `densify()`
infittisce il percorso (40 punti per segmento) per un'animazione fluida del trattore.
Passando sopra i focolai, il semaforo passa a **SPRUZZA** (dentro), **PRONTO** (<14 m),
**STOP** (lontano), calcolato da `nearestEdge()` (distanza dal bordo del focolaio più vicino).

---

## 4. NDVI e overlay del campo

Dall'ortomosaico multispettrale si ricavano **due immagini** sovrapponibili alla mappa.

`build_overlay()`:
1. **Lettura bande** a risoluzione ridotta (~1100 px lato lungo, per asset web leggeri).
   La banda **Rosso = banda 3**; il **NIR** si individua automaticamente come la banda
   più luminosa tra la 4 e la 5 (sul vigneto il NIR riflette molto) → qui **banda 5**.
2. **NDVI** = (NIR − Rosso) / (NIR + Rosso), pixel per pixel; si maschera il *nodata*.
3. **Overlay NDVI** (`ndvi.png`): NDVI mappato su una scala colore **RdYlGn**
   (rosso = vigore basso/suolo, verde = canopy sana).
4. **Overlay RGB true-color** (`field.png`): bande Rosso/Verde/Blu con **stretch ai
   percentili 2–98 %** + gamma, e **bilanciamento gray-world** per togliere la dominante
   di colore tipica delle camere multispettrali → sembra una foto aerea naturale.
5. **Allineamento**: NDVI, RGB e maschera vengono **riproiettati in EPSG:4326**; dai
   loro confini si ricavano i **`bounds`** geografici, che Leaflet usa per posizionare
   l'immagine esattamente sopra il campo (e quindi sotto i focolai).

→ output **`data/field.png`**, **`data/ndvi.png`** e i `bounds` in `meta.json`.
Nell'app il pulsante in alto a destra alterna **Drone RGB / NDVI / Satellite**.

---

## 5. Valutazioni e statistiche

Calcolate sui dati reali e mostrate nelle tre card della schermata Mappa:

- **Focolai Esca** = numero di cluster → **11**.
- **Area del vigneto** = pixel validi dell'ortomosaico × area pixel → **≈ 1,52 ha**.
- **Superficie da irrorare** = somma delle aree dei cerchi-focolaio → **≈ 1610 m²**.
- **Fitofarmaco risparmiato** = 1 − (area trattata / area campo) → **≈ 89,4 %**
  (rispetto a un trattamento "a tappeto" sull'intero vigneto).
- **NDVI medio**: campionando l'NDVI alla posizione di ogni vite →
  **malate 0,763 vs sane 0,769**: le viti malate hanno vigore leggermente inferiore,
  coerente con la fisiopatologia dell'Esca. È la "prova" che la diagnosi e il dato
  spettrale concordano.

→ tutto in **`data/meta.json`** (più centro campo, bounds, nome campo, attribuzione).

---

## 6. Come l'app usa i dati

`index.html`:
- carica i dati con `<script src="data/data.js">` (bundle `window.FITO_DATA`), così
  **funziona anche aperta col doppio clic** (`file://`); se servita via http usa il
  `fetch` come fallback;
- popola i `FOCI`, il `route`, le statistiche e le etichette (`loadRealData`/`applyRealUI`);
- centra le mappe sul vigneto reale e aggiunge gli overlay con `L.imageOverlay`;
- la sincronizzazione "Sincronizza Mappa Interventi" carica davvero questi dati.

---

## 7. Riprodurre la pipeline

```bash
pip install pyshp pyproj rasterio numpy Pillow requests
python tools/process_dataset.py
```

Lo script scarica i file da Zenodo in `raw/` (la prima volta ~3.2 GB) e rigenera tutto
`data/`. Parametri regolabili in cima allo script: `VINEYARD`, `LINK_DIST`,
`FOC_MARGIN`, `RAD_CAP`.

> `raw/` è in `.gitignore` (file grezzi pesanti, rigenerabili). Gli asset finali in
> `data/` **sono** versionati, così l'app funziona senza dover riscaricare nulla.

---

## 8. Limiti e note oneste

- I focolai derivano da una **diagnosi visiva a terra** (colonna `Esca`), non da una
  classificazione automatica dell'immagine: lo spettro NDVI è usato come **conferma**,
  non come rilevatore. Una versione futura potrebbe segmentare i focolai direttamente
  dall'NDVI/multispettrale.
- Il raggio dei focolai e la geometria del tragitto sono **stime ragionevoli** a fini
  dimostrativi, non un piano agronomico operativo.
- L'allineamento overlay↔focolai è sub-metrico ma non millimetrico (datum ETRS89 vs
  WGS84, riproiezione, decimazione dell'immagine).

---

## Attribuzione

Dati: **Vélez, Ariza-Sentís & Valente (2024)** — *Precision viticulture dataset*,
Zenodo [10.5281/zenodo.10362568](https://doi.org/10.5281/zenodo.10362568), **CC-BY-4.0**.

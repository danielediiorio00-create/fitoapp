#!/usr/bin/env python3
"""
FITOAPP — Pipeline di processing (una tantum).

Fonte: "Precision viticulture dataset ... Northern Spain, July 2022"
        Velez S., Ariza-Sentis M., Valente J.  CC-BY-4.0
        Zenodo: https://doi.org/10.5281/zenodo.10362568

Da `datasheet.csv` (una vite per riga, con coordinate UTM ETRS89/29N e
colonna `Esca` = YES/NO) ricava i PUNTI DI IRRORAZIONE REALI (= viti malate),
li raggruppa in focolai, costruisce un percorso a serpentina per la guida,
calcola un overlay NDVI reale dall'ortomosaico multispettrale e produce gli
asset web in ../data/.

Output:
  data/targets.geojson  focolai reali (centro lat/lng, raggio, n viti)
  data/route.geojson     percorso di guida (waypoint lat/lng)
  data/field.png         clip NDVI a colori dell'ortomosaico (overlay Leaflet)
  data/meta.json         bounds overlay, centro, statistiche, attribuzione
"""
import os, csv, json, math
from pyproj import Transformer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW  = os.path.join(ROOT, "raw")
DATA = os.path.join(ROOT, "data")
os.makedirs(DATA, exist_ok=True)

SRC_EPSG   = "EPSG:25829"          # ETRS89 / UTM 29N: datum delle coord. viti (dal .prj)
VINEYARD   = "B7"                  # il datasheet contiene B7 e B9; l'ortomosaico e' di B7
LINK_DIST  = 6.5                   # m: viti malate piu' vicine -> stesso focolaio
FOC_MARGIN = 2.0                   # m: margine raggio focolaio
RAD_CAP    = 14.0                  # m: raggio massimo focolaio (evita blob elongati)
ATTRIB = ("Vélez, Ariza-Sentís & Valente (2024) — Precision viticulture "
          "dataset, Zenodo 10362568, CC-BY-4.0")

to_wgs = Transformer.from_crs(SRC_EPSG, "EPSG:4326", always_xy=True)
def lonlat(x, y):
    lon, lat = to_wgs.transform(x, y)
    return lon, lat

# ── 1. Lettura viti + stato Esca ─────────────────────────────────────────────
def load_vines():
    path = os.path.join(RAW, "datasheet.csv")
    vines = []
    with open(path, encoding="utf-8-sig", errors="replace") as fh:
        for row in csv.DictReader(fh):
            if (row.get("Vineyard") or "").strip() != VINEYARD:
                continue
            try:
                x, y = float(row["X"]), float(row["Y"])
            except (TypeError, ValueError):
                continue
            esca = (row.get("Esca") or "").strip().upper()
            diseased = esca.startswith("YES")           # YES -> malata
            try:
                clusters = int(row.get("Number of grape clusters") or 0)
            except ValueError:
                clusters = 0
            vines.append({"x": x, "y": y, "diseased": diseased, "clusters": clusters})
    return vines

# ── 2. Clustering viti malate -> focolai (union-find su grafo di distanza) ────
def cluster(points, link):
    n = len(points)
    parent = list(range(n))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for i in range(n):
        for j in range(i + 1, n):
            dx = points[i]["x"] - points[j]["x"]
            dy = points[i]["y"] - points[j]["y"]
            if dx * dx + dy * dy <= link * link:
                union(i, j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(points[i])
    return list(groups.values())

def foci_from_clusters(groups):
    foci = []
    for g in groups:
        cx = sum(p["x"] for p in g) / len(g)
        cy = sum(p["y"] for p in g) / len(g)
        rad = max((math.hypot(p["x"] - cx, p["y"] - cy) for p in g), default=0.0)
        rad = round(min(rad + FOC_MARGIN, RAD_CAP), 1)
        lon, lat = lonlat(cx, cy)
        foci.append({"lat": round(lat, 7), "lng": round(lon, 7),
                     "r": rad, "n": len(g), "cx": cx, "cy": cy})
    foci.sort(key=lambda f: (-f["n"], f["lat"]))
    return foci

# ── 3. Geometria di campo (convex hull + area) ───────────────────────────────
def convex_hull(pts):
    pts = sorted(set((round(p["x"], 3), round(p["y"], 3)) for p in pts))
    if len(pts) < 3:
        return pts
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def poly_area(poly):
    s = 0.0
    for i in range(len(poly)):
        x1, y1 = poly[i]; x2, y2 = poly[(i + 1) % len(poly)]
        s += x1 * y2 - x2 * y1
    return abs(s) / 2.0

# ── 4. Percorso a serpentina lungo i filari (PCA sulle posizioni viti) ────────
def serpentine_route(vines, n_pass=5):
    cx = sum(v["x"] for v in vines) / len(vines)
    cy = sum(v["y"] for v in vines) / len(vines)
    sxx = sum((v["x"]-cx)**2 for v in vines)
    syy = sum((v["y"]-cy)**2 for v in vines)
    sxy = sum((v["x"]-cx)*(v["y"]-cy) for v in vines)
    # autovettore maggiore della matrice di covarianza 2x2 (direzione filari)
    theta = 0.5 * math.atan2(2 * sxy, sxx - syy)
    ux, uy = math.cos(theta), math.sin(theta)      # asse lungo (filare)
    vx, vy = -uy, ux                               # asse corto (tra filari)
    us = [ (v["x"]-cx)*ux + (v["y"]-cy)*uy for v in vines ]
    vs = [ (v["x"]-cx)*vx + (v["y"]-cy)*vy for v in vines ]
    umin, umax = min(us), max(us)
    vmin, vmax = min(vs), max(vs)
    pad = 4.0
    umin -= pad; umax += pad
    wp = []
    for i in range(n_pass):
        vv = vmin + (vmax - vmin) * (i / (n_pass - 1)) if n_pass > 1 else (vmin+vmax)/2
        ends = [umin, umax] if i % 2 == 0 else [umax, umin]
        for uu in ends:
            x = cx + uu * ux + vv * vx
            y = cy + uu * uy + vv * vy
            lon, lat = lonlat(x, y)
            wp.append({"lat": round(lat, 7), "lng": round(lon, 7)})
    return wp

# ── 5. Overlay NDVI dall'ortomosaico multispettrale ──────────────────────────
RDYLGN = [(0.00,(165,0,38)),(0.25,(244,109,67)),(0.50,(255,255,191)),
          (0.75,(166,217,106)),(1.00,(26,152,80))]
def ramp(t):
    t = max(0.0, min(1.0, t))
    for i in range(len(RDYLGN) - 1):
        t0, c0 = RDYLGN[i]; t1, c1 = RDYLGN[i + 1]
        if t <= t1:
            f = 0 if t1 == t0 else (t - t0) / (t1 - t0)
            return tuple(int(c0[k] + (c1[k] - c0[k]) * f) for k in range(3))
    return RDYLGN[-1][1]

def build_overlay(vines, max_dim=1100, ndvi_lo=0.10, ndvi_hi=0.85):
    """Ritorna (bounds_leaflet, ndvi_diseased_mean, ndvi_healthy_mean) o None."""
    import numpy as np
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import Affine, array_bounds
    from rasterio.warp import calculate_default_transform, reproject
    from PIL import Image

    tif = os.path.join(RAW, "ortho_b7.tif")
    if not os.path.exists(tif):
        print("  [overlay] ortho_b7.tif non presente — salto NDVI")
        return None
    try:
        return _build_overlay_impl(tif, vines, max_dim, ndvi_lo, ndvi_hi,
                                   np, rasterio, Resampling, Affine, array_bounds,
                                   calculate_default_transform, reproject, Image)
    except Exception as e:
        print(f"  [overlay] errore lettura ortomosaico ({e!r}) — salto NDVI")
        return None

def _build_overlay_impl(tif, vines, max_dim, ndvi_lo, ndvi_hi,
                        np, rasterio, Resampling, Affine, array_bounds,
                        calculate_default_transform, reproject, Image):
    with rasterio.open(tif) as src:
        print(f"  [overlay] bande={src.count} dtype={src.dtypes[0]} "
              f"size={src.width}x{src.height} nodata={src.nodata}")
        means = []
        for b in range(1, src.count + 1):
            s = src.read(b, out_shape=(256, 256), resampling=Resampling.average).astype("float64")
            means.append(float(s[s > 0].mean()) if (s > 0).any() else 0.0)
        print("  [overlay] media bande:", [round(m, 1) for m in means])
        # Red = banda 3 (B,G,R,...); NIR = banda piu' luminosa tra 4 e 5
        red_idx = 3 if src.count >= 3 else 1
        cand = [i for i in (4, 5) if i <= src.count]
        nir_idx = max(cand, key=lambda i: means[i - 1]) if cand else src.count
        print(f"  [overlay] red=banda{red_idx} nir=banda{nir_idx}")

        scale = max(src.width, src.height) / max_dim
        ow, oh = max(1, int(src.width / scale)), max(1, int(src.height / scale))
        red = src.read(red_idx, out_shape=(oh, ow), resampling=Resampling.average).astype("float32")
        nir = src.read(nir_idx, out_shape=(oh, ow), resampling=Resampling.average).astype("float32")
        grn = src.read(2, out_shape=(oh, ow), resampling=Resampling.average).astype("float32")
        blu = src.read(1, out_shape=(oh, ow), resampling=Resampling.average).astype("float32")
        t = src.transform * Affine.scale(src.width / ow, src.height / oh)
        valid = (red > 0) & (nir > 0)
        if src.nodata is not None:
            valid &= (red != src.nodata) & (nir != src.nodata)
        denom = nir + red
        ndvi = np.where(denom != 0, (nir - red) / denom, 0.0).astype("float32")

        # Area reale del blocco vigneto = pixel validi * area pixel (UTM, m)
        px_area = abs(t.a * t.e)
        field_area = float(valid.sum()) * px_area

        # NDVI medio viti malate vs sane (campionamento nel raster UTM decimato)
        def sample(v):
            col = int((v["x"] - t.c) / t.a); row = int((v["y"] - t.f) / t.e)
            if 0 <= row < oh and 0 <= col < ow and valid[row, col]:
                return float(ndvi[row, col])
            return None
        d_vals = [s for s in (sample(v) for v in vines if v["diseased"]) if s is not None]
        h_vals = [s for s in (sample(v) for v in vines if not v["diseased"]) if s is not None]
        ndvi_d = round(sum(d_vals) / len(d_vals), 3) if d_vals else None
        ndvi_h = round(sum(h_vals) / len(h_vals), 3) if h_vals else None

        # Reproiezione bande + mask -> EPSG:4326 (allineamento con i focolai)
        west, south, east, north = array_bounds(oh, ow, t)
        dst_crs = "EPSG:4326"
        dt, dw, dh = calculate_default_transform(src.crs, dst_crs, ow, oh,
                                                 west, south, east, north)

        def warp(arr, rs):
            out = np.zeros((dh, dw), "float32")
            reproject(arr.astype("float32"), out, src_transform=t, src_crs=src.crs,
                      dst_transform=dt, dst_crs=dst_crs, resampling=rs)
            return out
        R = warp(red, Resampling.bilinear)
        G = warp(grn, Resampling.bilinear)
        B = warp(blu, Resampling.bilinear)
        ndvi_d4 = warp(ndvi, Resampling.bilinear)
        mask_d4 = warp(valid, Resampling.nearest)

    h, w = ndvi_d4.shape
    alpha = np.where(mask_d4 > 0.5, 235, 0).astype("uint8")
    m = mask_d4 > 0.5

    # ── Overlay 1: true-color RGB con stretch percentile (foto aerea naturale) ──
    def stretch(band, gamma=1.15):
        v = band[m]
        lo, hi = (np.percentile(v, 2), np.percentile(v, 98)) if v.size else (0.0, 1.0)
        if hi <= lo:
            hi = lo + 1e-6
        s = np.clip((band - lo) / (hi - lo), 0, 1) ** (1.0 / gamma)
        return (s * 255).astype("uint8")
    chans = [stretch(R).astype("float32"), stretch(G).astype("float32"), stretch(B).astype("float32")]
    # gray-world: riallinea le medie dei canali per togliere la dominante di colore
    means_c = [c[m].mean() if m.any() else 1.0 for c in chans]
    gray = sum(means_c) / 3.0
    chans = [np.clip(c * (gray / mc), 0, 255) if mc > 0 else c for c, mc in zip(chans, means_c)]
    rgb = np.zeros((h, w, 4), "uint8")
    for k in range(3):
        rgb[..., k] = chans[k].astype("uint8")
    rgb[..., 3] = alpha
    Image.fromarray(rgb, "RGBA").save(os.path.join(DATA, "field.png"), optimize=True)

    # ── Overlay 2: mappa NDVI a colori (RdYlGn) ──
    nd = np.zeros((h, w, 4), "uint8")
    tnorm = np.clip((ndvi_d4 - ndvi_lo) / (ndvi_hi - ndvi_lo), 0, 1)
    lut = np.array([ramp(i / 255.0) for i in range(256)], "uint8")
    idx = (tnorm * 255).astype("uint8")
    nd[..., 0] = lut[idx, 0]; nd[..., 1] = lut[idx, 1]; nd[..., 2] = lut[idx, 2]
    nd[..., 3] = alpha
    Image.fromarray(nd, "RGBA").save(os.path.join(DATA, "ndvi.png"), optimize=True)
    print(f"  [overlay] field.png + ndvi.png salvati {w}x{h}")

    wN, sN, eN, nN = array_bounds(h, w, dt)
    bounds = [[round(sN, 7), round(wN, 7)], [round(nN, 7), round(eN, 7)]]
    return bounds, ndvi_d, ndvi_h, round(field_area)

# ── main ─────────────────────────────────────────────────────────────────────
def main():
    vines = load_vines()
    diseased = [v for v in vines if v["diseased"]]
    print(f"Viti totali: {len(vines)} | malate (Esca=YES): {len(diseased)}")

    foci = foci_from_clusters(cluster(diseased, LINK_DIST))
    print(f"Focolai (cluster a {LINK_DIST} m): {len(foci)}")

    treated_area = sum(math.pi * f["r"] ** 2 for f in foci)

    cx = sum(v["x"] for v in vines) / len(vines)
    cy = sum(v["y"] for v in vines) / len(vines)
    clon, clat = lonlat(cx, cy)

    overlay = build_overlay(vines)
    bounds = overlay[0] if overlay else None
    ndvi_d = overlay[1] if overlay else None
    ndvi_h = overlay[2] if overlay else None
    # area campo: estensione reale dell'ortomosaico se disponibile, altrimenti hull viti
    field_area = overlay[3] if overlay else poly_area(convex_hull(vines))
    saved_pct = round((1 - treated_area / field_area) * 100, 1) if field_area else 0.0

    # targets.geojson
    tg = {"type": "FeatureCollection",
          "features": [{"type": "Feature",
                        "geometry": {"type": "Point", "coordinates": [f["lng"], f["lat"]]},
                        "properties": {"r": f["r"], "vines": f["n"]}} for f in foci]}
    json.dump(tg, open(os.path.join(DATA, "targets.geojson"), "w"), indent=1)

    # route.geojson
    wp = serpentine_route(vines)
    rt = {"type": "FeatureCollection",
          "features": [{"type": "Feature",
                        "geometry": {"type": "LineString",
                                     "coordinates": [[p["lng"], p["lat"]] for p in wp]},
                        "properties": {"role": "spray-path"}}]}
    json.dump(rt, open(os.path.join(DATA, "route.geojson"), "w"), indent=1)

    meta = {
        "field": "Vigneto B7 — Galizia, Spagna",
        "center": {"lat": round(clat, 7), "lng": round(clon, 7)},
        "bounds": bounds,
        "stats": {
            "vines_total": len(vines),
            "vines_diseased": len(diseased),
            "foci": len(foci),
            "field_area_m2": round(field_area),
            "treated_area_m2": round(treated_area),
            "saved_pct": saved_pct,
            "ndvi_diseased": ndvi_d,
            "ndvi_healthy": ndvi_h,
        },
        "source": {"name": "Precision viticulture dataset — Northern Spain (2022)",
                   "doi": "10.5281/zenodo.10362568", "license": "CC-BY-4.0",
                   "attribution": ATTRIB},
    }
    json.dump(meta, open(os.path.join(DATA, "meta.json"), "w"), indent=1)

    # Bundle JS: permette all'app di funzionare anche aperta via file:// (doppio
    # clic), dove il fetch dei file locali e' bloccato dal browser.
    with open(os.path.join(DATA, "data.js"), "w", encoding="utf-8") as fh:
        fh.write("window.FITO_DATA=" +
                 json.dumps({"meta": meta, "targets": tg, "route": rt}, ensure_ascii=False) +
                 ";\n")

    print("\n=== RISULTATO ===")
    print(json.dumps(meta["stats"], indent=1, ensure_ascii=False))
    print("Output scritto in", DATA)

if __name__ == "__main__":
    main()

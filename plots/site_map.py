"""Redraw the study-area map from the cached regional layers and the S1 station manifest."""
from __future__ import annotations
import argparse
import gzip
import json
import math
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LightSource
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle, Polygon as MplPolygon
import numpy as np
import pandas as pd
from PIL import Image
from pyproj import Geod
from shapely import make_valid
from shapely.geometry import MultiPolygon, Polygon, shape
from shapely.ops import unary_union

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--layers",type=Path,required=True)
parser.add_argument("--s1",type=Path,required=True)
parser.add_argument("--out",type=Path,required=True)
args=parser.parse_args()
DATA,OUTPUTS=args.layers,args.out
settings=json.loads((DATA/"map_settings.json").read_text())
ZBY_BBOX,DXC_BBOX,REGIONAL_BBOX=(tuple(settings[k]) for k in ("zby_bbox","dxc_bbox","regional_bbox"))
ZBY,DXC=(tuple(settings[k]) for k in ("zby_locator","dxc_locator"))
GEOD=Geod(ellps="WGS84")

LAKE_FACE = "#8EC9D8"
LAKE_EDGE = "#155E75"
STREAM_COLOR = "#438EB2"
GRASS_COLOR = "#789B4A"
CONTOUR_COLOR = "#61594D"
PRIMARY_COLOR = "#B32635"
COMPARISON_COLOR = "#C47C2C"

def read_single_tile_dem(tile: str, bbox: tuple[float, float, float, float]):
    west, south, east, north = bbox
    lat0 = int(tile[1:3])
    lon0 = int(tile[4:7])
    if tile[0] == "S":
        lat0 *= -1
    if tile[3] == "W":
        lon0 *= -1
    array = np.frombuffer(gzip.decompress((DATA / f"{tile}.hgt.gz").read_bytes()), dtype=">i2").reshape(3601, 3601)
    row0 = max(0, int(math.floor((lat0 + 1 - north) * 3600)))
    row1 = min(3601, int(math.ceil((lat0 + 1 - south) * 3600)) + 1)
    col0 = max(0, int(math.floor((west - lon0) * 3600)))
    col1 = min(3601, int(math.ceil((east - lon0) * 3600)) + 1)
    dem = array[row0:row1, col0:col1].astype(float)
    dem[dem <= -32768] = np.nan
    if np.isnan(dem).any():
        dem[np.isnan(dem)] = np.nanmedian(dem)
    extent = (
        lon0 + col0 / 3600,
        lat0 + 1 - (row1 - 1) / 3600,
        lon0 + (col1 - 1) / 3600,
        lat0 + 1 - row0 / 3600,
    )
    return dem, extent

def read_tiled_dem(
    tiles: tuple[str, ...],
    bbox: tuple[float, float, float, float],
):
    """Read adjacent one-degree SRTM tiles and remove their duplicated edge columns."""
    parts = [read_single_tile_dem(tile, bbox) for tile in tiles]
    parts.sort(key=lambda item: item[1][0])
    row_count = parts[0][0].shape[0]
    south = parts[0][1][1]
    north = parts[0][1][3]
    for array, extent in parts[1:]:
        if array.shape[0] != row_count or not np.isclose(extent[1], south) or not np.isclose(extent[3], north):
            raise RuntimeError("Adjacent SRTM tiles do not share the same latitude grid")
    arrays = [parts[0][0]]
    for array, _extent in parts[1:]:
        arrays.append(array[:, 1:])
    mosaic = np.hstack(arrays)
    extent = (parts[0][1][0], south, parts[-1][1][2], north)
    return mosaic, extent

def tile_xy(lon: float, lat: float, zoom: int) -> tuple[int, int]:
    n = 2**zoom
    x = int((lon + 180.0) / 360.0 * n)
    lat_rad = math.radians(lat)
    y = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
    return x, y

def tile_lat(y: float, zoom: int) -> float:
    n = 2**zoom
    return math.degrees(math.atan(math.sinh(math.pi * (1.0 - 2.0 * y / n))))

def get_regional_relief() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    zoom = 5
    west, south, east, north = REGIONAL_BBOX
    x0, y0 = tile_xy(west, north, zoom)
    x1, y1 = tile_xy(east, south, zoom)
    tiles: list[list[np.ndarray]] = []
    for y in range(y0, y1 + 1):
        row = []
        for x in range(x0, x1 + 1):
            path = DATA / f"terrarium_z{zoom}_{x}_{y}.png"
            rgb = np.asarray(Image.open(path).convert("RGB"), dtype=float)
            elevation = rgb[:, :, 0] * 256.0 + rgb[:, :, 1] + rgb[:, :, 2] / 256.0 - 32768.0
            row.append(elevation)
        tiles.append(row)
    relief = np.block(tiles)
    height, width = relief.shape
    n = 2**zoom
    lon_edges = np.linspace(x0 / n * 360.0 - 180.0, (x1 + 1) / n * 360.0 - 180.0, width + 1)
    y_edges = np.linspace(y0, y1 + 1, height + 1)
    lat_edges = np.array([tile_lat(y, zoom) for y in y_edges])
    return relief, lon_edges, lat_edges

def add_polygon(ax, polygon: Polygon, **kwargs) -> None:
    ax.add_patch(MplPolygon(np.asarray(polygon.exterior.coords), closed=True, **kwargs))
    for ring in polygon.interiors:
        ax.add_patch(
            MplPolygon(
                np.asarray(ring.coords),
                closed=True,
                facecolor=ax.get_facecolor(),
                edgecolor="none",
                zorder=kwargs.get("zorder", 1) + 0.1,
            )
        )

def iter_polygons(geometry):
    if isinstance(geometry, Polygon):
        yield geometry
    elif isinstance(geometry, MultiPolygon):
        yield from geometry.geoms
    elif hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            yield from iter_polygons(part)

def add_geometry_exterior_outline(ax, geometry, *, color: str, linewidth: float, zorder: int) -> None:
    """Draw only exterior rings so a transparent national outline cannot mask province fills."""
    for polygon in iter_polygons(geometry):
        exterior = np.asarray(polygon.exterior.coords)
        ax.plot(
            exterior[:, 0],
            exterior[:, 1],
            color=color,
            linewidth=linewidth,
            solid_joinstyle="round",
            solid_capstyle="round",
            zorder=zorder,
        )

def format_lon(value: float, _position=None) -> str:
    return f"{value:.1f}°E"

def format_lat(value: float, _position=None) -> str:
    return f"{value:.1f}°N"

def add_north_arrow(ax, x: float = 0.94, y: float = 0.90, size: float = 0.12) -> None:
    ax.annotate(
        "N",
        xy=(x, y),
        xytext=(x, y - size),
        xycoords="axes fraction",
        textcoords="axes fraction",
        ha="center",
        va="center",
        fontsize=7.2,
        fontweight="bold",
        color="#202020",
        arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#202020", shrinkA=3, shrinkB=0),
        zorder=30,
    )

def add_scale_bar(
    ax,
    lon: float,
    lat: float,
    length_km: float,
    bbox: tuple[float, float, float, float],
    segments: int = 2,
) -> None:
    segment_km = length_km / segments
    current_lon = lon
    bar_height = (bbox[3] - bbox[1]) * 0.018
    for index in range(segments):
        next_lon, _, _ = GEOD.fwd(current_lon, lat, 90, segment_km * 1000)
        ax.add_patch(
            Rectangle(
                (current_lon, lat),
                next_lon - current_lon,
                bar_height,
                facecolor="#202020" if index % 2 == 0 else "white",
                edgecolor="#202020",
                linewidth=0.6,
                zorder=30,
            )
        )
        current_lon = next_lon
    label_y = lat - bar_height * 1.3
    ax.text(lon, label_y, "0", ha="center", va="top", fontsize=5.5, zorder=31)
    ax.text(current_lon, label_y, f"{length_km:g} km", ha="center", va="top", fontsize=5.5, zorder=31)

def sample_dem(dem: np.ndarray, extent, lon: float, lat: float) -> float:
    west, south, east, north = extent
    col = int(round((lon - west) / (east - west) * (dem.shape[1] - 1)))
    row = int(round((north - lat) / (north - south) * (dem.shape[0] - 1)))
    row = max(0, min(dem.shape[0] - 1, row))
    col = max(0, min(dem.shape[1] - 1, col))
    return float(dem[row, col])

def add_terrain_flow_arrow(ax, visible: np.ndarray, dem: np.ndarray, extent) -> None:
    if len(visible) < 4:
        return
    adjusted = visible.copy()
    adjusted[:, 0] *= math.cos(math.radians(float(np.mean(visible[:, 1]))))
    lengths = np.sqrt(np.sum(np.diff(adjusted, axis=0) ** 2, axis=1))
    cumulative = np.concatenate(([0.0], np.cumsum(lengths)))
    total = float(cumulative[-1])
    if total < 0.035:
        return
    start_probe = visible[max(0, int(len(visible) * 0.08))]
    end_probe = visible[min(len(visible) - 1, int(len(visible) * 0.92))]
    if sample_dem(dem, extent, *start_probe) < sample_dem(dem, extent, *end_probe):
        visible = visible[::-1].copy()
        adjusted = visible.copy()
        adjusted[:, 0] *= math.cos(math.radians(float(np.mean(visible[:, 1]))))
        lengths = np.sqrt(np.sum(np.diff(adjusted, axis=0) ** 2, axis=1))
        cumulative = np.concatenate(([0.0], np.cumsum(lengths)))
        total = float(cumulative[-1])

    def point_along(distance: float) -> np.ndarray:
        index = int(np.searchsorted(cumulative, distance, side="right") - 1)
        index = max(0, min(index, len(visible) - 2))
        span = cumulative[index + 1] - cumulative[index]
        fraction = 0 if span <= 0 else (distance - cumulative[index]) / span
        return visible[index] + fraction * (visible[index + 1] - visible[index])

    center = total * 0.58
    tail = point_along(max(0, center - min(0.012, total * 0.08)))
    head = point_along(min(total, center + min(0.012, total * 0.08)))
    ax.annotate(
        "",
        xy=head,
        xytext=tail,
        arrowprops=dict(
            arrowstyle="-|>",
            color="#287FA8",
            lw=0.7,
            mutation_scale=6,
            shrinkA=0,
            shrinkB=0,
        ),
        zorder=8,
    )

def plot_streams(ax, payload: dict, bbox, dem, extent) -> None:
    west, south, east, north = bbox
    for element in payload.get("elements", []):
        geometry = element.get("geometry", [])
        if len(geometry) < 2:
            continue
        coords = np.array([(node["lon"], node["lat"]) for node in geometry], dtype=float)
        inside = (
            (coords[:, 0] >= west - 0.01)
            & (coords[:, 0] <= east + 0.01)
            & (coords[:, 1] >= south - 0.01)
            & (coords[:, 1] <= north + 0.01)
        )
        indices = np.flatnonzero(inside)
        if len(indices) < 2:
            continue
        visible = coords[indices[0] : indices[-1] + 1]
        kind = element.get("tags", {}).get("waterway", "stream")
        is_primary = kind == "river" or len(visible) > 70
        ax.plot(
            visible[:, 0],
            visible[:, 1],
            color=STREAM_COLOR,
            linewidth=0.72 if is_primary else 0.36,
            alpha=0.82 if is_primary else 0.60,
            zorder=6,
        )
        if is_primary:
            add_terrain_flow_arrow(ax, visible, dem, extent)

def grassland_rgba(path: Path) -> np.ndarray:
    rgb = np.asarray(Image.open(path).convert("RGB"), dtype=float) / 255.0
    rgba = np.zeros((*rgb.shape[:2], 4), dtype=float)
    rgba[:, :, :3] = np.array(mpl.colors.to_rgb(GRASS_COLOR))
    rgba[:, :, 3] = np.where(rgb.max(axis=2) > 0.18, 0.42, 0.0)
    return rgba

def plot_lake_context(
    ax,
    *,
    dem: np.ndarray,
    dem_extent,
    bbox,
    lake_payload: dict,
    streams_payload: dict,
    grassland_path: Path,
    contour_interval: int,
    title: str,
    panel_label: str,
    scale_km: float,
    scale_origin: tuple[float, float],
):
    terrain_cmap = LinearSegmentedColormap.from_list(
        "plateau_terrain",
        ["#EEE8DC", "#D7CCB9", "#BDAE96", "#A08D75", "#E5E0D7", "#F7F7F4"],
    )
    west, south, east, north = dem_extent
    mean_lat = (south + north) / 2
    lon_res = (east - west) / (dem.shape[1] - 1)
    lat_res = (north - south) / (dem.shape[0] - 1)
    dx = lon_res * 111320 * math.cos(math.radians(mean_lat))
    dy = lat_res * 110574
    hillshade = LightSource(azdeg=315, altdeg=42).hillshade(dem, vert_exag=1.0, dx=dx, dy=dy)
    raster = ax.imshow(
        dem,
        extent=(west, east, south, north),
        origin="upper",
        cmap=terrain_cmap,
        vmin=4380,
        vmax=5600,
        interpolation="bilinear",
        rasterized=True,
        zorder=1,
    )
    ax.imshow(
        hillshade,
        extent=(west, east, south, north),
        origin="upper",
        cmap="gray",
        vmin=0,
        vmax=1,
        alpha=0.25,
        interpolation="bilinear",
        rasterized=True,
        zorder=2,
    )
    ax.imshow(
        grassland_rgba(grassland_path),
        extent=(bbox[0],bbox[2],bbox[1],bbox[3]),
        origin="upper",
        interpolation="nearest",
        rasterized=True,
        zorder=2.5,
    )
    lon_values = np.linspace(west, east, dem.shape[1])
    lat_values = np.linspace(north, south, dem.shape[0])
    levels = np.arange(
        math.ceil(np.nanmin(dem) / contour_interval) * contour_interval,
        math.floor(np.nanmax(dem) / contour_interval) * contour_interval + contour_interval,
        contour_interval,
    )
    contours = ax.contour(
        lon_values,
        lat_values,
        dem,
        levels=levels,
        colors=CONTOUR_COLOR,
        linewidths=0.28,
        alpha=0.46,
        zorder=3,
    )
    if len(contours.levels):
        ax.clabel(contours, levels=contours.levels[::2], inline=True, fontsize=4.0, fmt="%d")
    plot_streams(ax, streams_payload, bbox, dem, dem_extent)
    lake_geometry = shape(lake_payload["features"][0]["geometry"])
    for polygon in iter_polygons(lake_geometry):
        add_polygon(
            ax,
            polygon,
            facecolor=LAKE_FACE,
            edgecolor=LAKE_EDGE,
            linewidth=0.85,
            alpha=0.97,
            zorder=9,
        )
    ax.set_xlim(bbox[0], bbox[2])
    ax.set_ylim(bbox[1], bbox[3])
    x_step = 0.1 if (bbox[2] - bbox[0]) > 0.35 else 0.05
    y_step = 0.1 if bbox == ZBY_BBOX else 0.05
    ax.set_xticks(np.arange(math.ceil(bbox[0] / x_step) * x_step, bbox[2] + 0.001, x_step))
    ax.set_yticks(np.arange(math.ceil(bbox[1] / y_step) * y_step, bbox[3] + 0.001, y_step))
    ax.xaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(
            lambda value, position: f"{value:.1f}°E"
        )
    )
    ax.yaxis.set_major_formatter(
        mpl.ticker.FuncFormatter(
            lambda value, position: f"{value:.1f}°N"
            if bbox == ZBY_BBOX
            else f"{value:.2f}°N"
        )
    )
    ax.grid(color="white", linewidth=0.40, alpha=0.60, zorder=4)
    ax.tick_params(labelsize=5.5, length=2.0, pad=1.2)
    ax.set_aspect(1 / math.cos(math.radians(mean_lat)))
    ax.text(
        0.012,
        0.986,
        panel_label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="none", pad=1.1, alpha=0.88),
        zorder=40,
    )
    ax.text(
        0.065,
        0.975,
        title,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.9,
        fontweight="bold",
        color="#173C49",
        bbox=dict(
            boxstyle="round,pad=0.18",
            facecolor="white",
            edgecolor="#B8B8B8",
            linewidth=0.35,
            alpha=0.88,
        ),
        zorder=39,
    )
    add_north_arrow(ax)
    add_scale_bar(ax, *scale_origin, scale_km, bbox)
    return raster

def build_figure() -> list[Path]:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    def layer(filename):
        return json.loads((DATA/filename).read_text(encoding="utf8"))
    china_provinces=layer("china_provinces_adcode_100000.geojson")
    zby_dem,zby_extent=read_tiled_dem(("N31E083","N31E084"), ZBY_BBOX)
    dxc_dem,dxc_extent=read_tiled_dem(("N31E086","N31E087"), DXC_BBOX)
    zby_lake=layer("zabuye_lake_osm_relation_11763855.geojson")
    dxc_lake=layer("dxc_lake_osm_relation_2309319.geojson")
    zby_streams=layer("zabuye_waterways_osm_overpass.json")
    dxc_streams=layer("dxc_waterways_osm_overpass_v3.json")
    regional_dem,regional_lon_edges,regional_lat_edges=get_regional_relief()
    zby_grass = DATA / "zabuye_worldcover_grassland_2021.png"
    dxc_grass = DATA / "dxc_worldcover_grassland_2021.png"
    for path in (zby_grass, dxc_grass):
        if not path.exists():
            raise RuntimeError(f"Missing Earth Engine layer: {path}")

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 7.2,
            "axes.linewidth": 0.68,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(183 / 25.4, 112 / 25.4), facecolor="white")
    fig.subplots_adjust(left=0.035, right=0.985, bottom=0.085, top=0.97)
    grid = fig.add_gridspec(
        6,
        2,
        width_ratios=(0.60, 1.62),
        height_ratios=(1, 1, 1, 1, 1, 1),
        wspace=0.065,
        hspace=0.22,
    )
    ax_china = fig.add_subplot(grid[0:3, 0])
    ax_region = fig.add_subplot(grid[3:6, 0])
    right_grid = grid[:, 1].subgridspec(
        2,
        1,
        height_ratios=(ZBY_BBOX[3] - ZBY_BBOX[1], DXC_BBOX[3] - DXC_BBOX[1]),
        hspace=0.10,
    )
    ax_zby = fig.add_subplot(right_grid[0, 0])
    ax_dxc = fig.add_subplot(right_grid[1, 0])

    province_geometries = []
    for feature in china_provinces["features"]:
        geometry = make_valid(shape(feature["geometry"]))
        province_geometries.append(geometry)
        is_tibet = str(feature.get("properties", {}).get("adcode", "")) == "540000"
        for polygon in iter_polygons(geometry):
            add_polygon(
                ax_china,
                polygon,
                facecolor="#D95B59" if is_tibet else "#D2D0CA",
                edgecolor="#A52A2A" if is_tibet else "#77736C",
                linewidth=0.70 if is_tibet else 0.40,
                alpha=0.96,
                zorder=4 if not is_tibet else 7,
            )
    national_geometry = unary_union(province_geometries)
    add_geometry_exterior_outline(
        ax_china,
        national_geometry,
        color="#343434",
        linewidth=0.68,
        zorder=8,
    )
    south_china_sea = ax_china.inset_axes([0.86, -0.08, 0.22, 0.33])
    south_china_sea.set_facecolor("white")
    for feature in china_provinces["features"]:
        geometry = make_valid(shape(feature["geometry"]))
        for polygon in iter_polygons(geometry):
            add_polygon(
                south_china_sea,
                polygon,
                facecolor="#D2D0CA",
                edgecolor="#77736C",
                linewidth=0.34,
                alpha=0.98,
                zorder=4,
            )
    add_geometry_exterior_outline(
        south_china_sea,
        national_geometry,
        color="#343434",
        linewidth=0.52,
        zorder=8,
    )
    south_china_sea.set_xlim(105.0, 125.0)
    south_china_sea.set_ylim(3.0, 25.5)
    south_china_sea.set_aspect(1 / math.cos(math.radians(15.0)))
    south_china_sea.set_xticks([])
    south_china_sea.set_yticks([])
    for spine in south_china_sea.spines.values():
        spine.set_linewidth(0.45)
        spine.set_color("#777777")
    ax_china.set_xlim(73.0, 135.5)
    ax_china.set_ylim(17.0, 54.5)
    ax_china.set_aspect(1 / math.cos(math.radians(35.0)))
    ax_china.text(
        0.50,
        0.985,
        "China",
        transform=ax_china.transAxes,
        ha="center",
        va="top",
        fontsize=7.2,
        fontweight="bold",
        color="#2D2D2D",
        zorder=12,
    )
    ax_china.text(
        88.5,
        31.4,
        "Tibet Autonomous\nRegion",
        ha="center",
        va="center",
        fontsize=5.2,
        fontweight="bold",
        color="white",
        linespacing=0.95,
        path_effects=[path_effects.withStroke(linewidth=1.0, foreground="#9E2828")],
        zorder=15,
    )
    ax_china.add_patch(
        Rectangle(
            (REGIONAL_BBOX[0], REGIONAL_BBOX[1]),
            REGIONAL_BBOX[2] - REGIONAL_BBOX[0],
            REGIONAL_BBOX[3] - REGIONAL_BBOX[1],
            facecolor="none",
            edgecolor="#111111",
            linewidth=0.80,
            joinstyle="miter",
            zorder=14,
        )
    )
    ax_china.text(
        0.02,
        0.03,
        settings["standard_map_review"],
        transform=ax_china.transAxes,
        fontsize=4.8,
        ha="left",
        va="bottom",
        bbox=dict(
            boxstyle="round,pad=0.16",
            facecolor="white",
            edgecolor="#B5B5B5",
            linewidth=0.35,
            alpha=0.90,
        ),
    )
    ax_china.set_axis_off()
    ax_china.set_anchor("S")
    ax_china.text(
        0.01,
        0.99,
        "a",
        transform=ax_china.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="none", pad=1.1, alpha=0.90),
    )

    regional_cmap = LinearSegmentedColormap.from_list(
        "regional_relief",
        ["#D9E1CF", "#D8CEB7", "#B7AA8E", "#978872", "#E7E4DE", "#FFFFFF"],
    )
    regional_hs = LightSource(azdeg=315, altdeg=45).hillshade(
        regional_dem, vert_exag=1.0, dx=4000, dy=4000
    )
    regional_raster = ax_region.pcolormesh(
        regional_lon_edges,
        regional_lat_edges,
        regional_dem,
        cmap=regional_cmap,
        vmin=0,
        vmax=6500,
        shading="auto",
        rasterized=True,
    )
    ax_region.pcolormesh(
        regional_lon_edges,
        regional_lat_edges,
        regional_hs,
        cmap="gray",
        vmin=0,
        vmax=1,
        alpha=0.24,
        shading="auto",
        rasterized=True,
    )
    ax_region.text(
        0.50,
        0.10,
        "Tibetan Plateau",
        transform=ax_region.transAxes,
        ha="center",
        va="center",
        fontsize=7.2,
        color="#5F5B52",
        alpha=0.80,
    )
    ax_region.scatter(
        *ZBY,
        marker="*",
        s=46,
        c=PRIMARY_COLOR,
        edgecolors="white",
        linewidths=0.55,
        zorder=12,
    )
    ax_region.scatter(
        *DXC,
        marker="D",
        s=22,
        c=COMPARISON_COLOR,
        edgecolors="white",
        linewidths=0.55,
        zorder=12,
    )
    ax_region.text(
        ZBY[0] - 0.25,
        ZBY[1] - 0.42,
        "ZBY",
        fontsize=6.1,
        fontweight="bold",
        color="#76131F",
        ha="right",
        va="top",
        zorder=13,
    )
    ax_region.text(
        DXC[0] + 0.30,
        DXC[1] + 0.22,
        "DXC",
        fontsize=6.1,
        fontweight="bold",
        color="#7E4A15",
        ha="left",
        va="bottom",
        zorder=13,
    )
    for bbox, color in ((ZBY_BBOX, PRIMARY_COLOR), (DXC_BBOX, COMPARISON_COLOR)):
        ax_region.add_patch(
            Rectangle(
                (bbox[0], bbox[1]),
                bbox[2] - bbox[0],
                bbox[3] - bbox[1],
                facecolor="none",
                edgecolor=color,
                linewidth=0.8,
                zorder=11,
            )
        )
    ax_region.set_xlim(REGIONAL_BBOX[0], REGIONAL_BBOX[2])
    ax_region.set_ylim(REGIONAL_BBOX[1], REGIONAL_BBOX[3])
    ax_region.set_xticks(np.arange(78, 93, 3))
    ax_region.set_yticks(np.arange(28, 37, 2))
    ax_region.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, p: f"{x:.0f}°E"))
    ax_region.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, p: f"{y:.0f}°N"))
    ax_region.grid(color="white", linewidth=0.35, alpha=0.64)
    ax_region.tick_params(labelsize=5.3, length=2, pad=1.1)
    ax_region.set_aspect(1 / math.cos(math.radians(31.8)))
    ax_region.set_anchor("N")
    ax_region.text(
        0.012,
        0.986,
        "b",
        transform=ax_region.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        fontweight="bold",
        bbox=dict(facecolor="white", edgecolor="none", pad=1.1, alpha=0.88),
        zorder=30,
    )
    add_north_arrow(ax_region, x=0.91, y=0.89, size=0.16)
    add_scale_bar(
        ax_region,
        85.5,
        27.85,
        500,
        REGIONAL_BBOX,
        segments=2,
    )
    plot_lake_context(
        ax_zby,
        dem=zby_dem,
        dem_extent=zby_extent,
        bbox=ZBY_BBOX,
        lake_payload=zby_lake,
        streams_payload=zby_streams,
        grassland_path=zby_grass,
        contour_interval=200,
        title="Zabuye Salt Lake (ZBY) · primary study lake",
        panel_label="c",
        scale_km=10,
        scale_origin=(83.825, 31.198),
    )
    ax_zby.text(
        84.050,
        31.457,
        "North Lake",
        fontsize=7.2,
        fontweight="bold",
        color="#0B4256",
        ha="center",
        va="center",
        zorder=20,
    )
    ax_zby.text(
        84.020,
        31.325,
        "South Lake",
        fontsize=7.2,
        fontweight="bold",
        color="#0B4256",
        ha="center",
        va="center",
        zorder=20,
    )
    river_label = ax_zby.text(
        84.245,
        31.528,
        "Luoju Zangbo",
        fontsize=5.8,
        fontstyle="italic",
        color="#1D6689",
        rotation=-31,
        rotation_mode="anchor",
        ha="center",
        va="bottom",
        zorder=20,
    )
    river_label.set_path_effects(
        [path_effects.withStroke(linewidth=1.7, foreground="white", alpha=0.88)]
    )

    plot_lake_context(
        ax_dxc,
        dem=dxc_dem,
        dem_extent=dxc_extent,
        bbox=DXC_BBOX,
        lake_payload=dxc_lake,
        streams_payload=dxc_streams,
        grassland_path=dxc_grass,
        contour_interval=200,
        title="Dangxiong Co (DXC) · comparison lake",
        panel_label="d",
        scale_km=10,
        scale_origin=(86.500, 31.472),
    )
    ax_dxc.text(
        86.748,
        31.569,
        "Dangxiong Co",
        fontsize=6.8,
        fontweight="bold",
        color="#0B4256",
        ha="center",
        va="center",
        zorder=20,
    )

    fig.add_artist(
        ConnectionPatch(
            xyA=(0.98, 0.72),
            coordsA=ax_region.transAxes,
            xyB=(0.0, 0.22),
            coordsB=ax_zby.transAxes,
            arrowstyle="-",
            color=PRIMARY_COLOR,
            linewidth=0.55,
            alpha=0.75,
            zorder=45,
        )
    )
    fig.add_artist(
        ConnectionPatch(
            xyA=(0.98, 0.56),
            coordsA=ax_region.transAxes,
            xyB=(0.0, 0.52),
            coordsB=ax_dxc.transAxes,
            arrowstyle="-",
            color=COMPARISON_COLOR,
            linewidth=0.55,
            alpha=0.75,
            zorder=45,
        )
    )

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=LAKE_FACE, edgecolor=LAKE_EDGE, linewidth=0.8),
        Rectangle((0, 0), 1, 1, facecolor=GRASS_COLOR, edgecolor="none", alpha=0.62),
        Line2D([0], [0], color=STREAM_COLOR, linewidth=0.8, marker=">", markersize=3.2),
        Line2D([0], [0], color=CONTOUR_COLOR, linewidth=0.4),
    ]
    fig.legend(
        legend_handles,
        ["Lake", "Grassland", "Watercourse and terrain-derived flow", "Elevation contour"],
        loc="lower center",
        bbox_to_anchor=(0.66, 0.001),
        ncol=4,
        frameon=False,
        fontsize=5.2,
        handlelength=1.5,
        columnspacing=1.1,
    )
    fig.text(
        0.995,
        0.006,
        "China locator: thematic rendering referenced to GS(2023)2766 | Relief: NASA SRTMGL1 | Land cover: ESA WorldCover 2021 | Lake and watercourses: © OpenStreetMap contributors",
        ha="right",
        va="bottom",
        fontsize=4.45,
        color="#4A4A4A",
    )

    points=pd.read_excel(args.s1,sheet_name="Sample_manifest_35")
    points=points[(points.matrix=="water") & (points.pooling_status=="Not pooled")]
    for lake,axis in (("ZBY",ax_zby),("DXC",ax_dxc)):
        for station,rows in points[points.lake==lake].groupby("resolved_station"):
            coords=rows[["longitude_degE","latitude_degN"]].drop_duplicates().dropna()
            if len(coords)!=1:raise ValueError(f"Ambiguous station coordinates for {lake} {station}")
            lon,lat=coords.iloc[0]
            label=station + ("*" if lake=="ZBY" and set(rows.season_or_campaign.str.lower())=={"autumn"} else "")
            axis.scatter(lon,lat,c="#a51b2f",edgecolors="white",linewidths=.3,s=8,zorder=50)
            axis.annotate(label,(lon,lat),xytext=(2,-4),textcoords="offset points",fontsize=4.7,zorder=51)
    paths=[]
    for suffix in ("png","svg","pdf"):
        path=OUTPUTS/("Figure_1a_map."+suffix)
        fig.savefig(path,dpi=300,bbox_inches="tight",pad_inches=.035)
        paths.append(path)
    plt.close(fig)
    return paths


if __name__=="__main__":
    for path in build_figure(): print(path)

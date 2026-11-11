"""Render the six cube-face views of a protein into one PNG.

The protein is drawn with cellscape's per-residue outlines, coloured by chain,
from the front, left, right, back, top and bottom, laid out in a row at a
common scale.

"Front" is chosen automatically by aligning the protein's principal axes with
the screen: its longest dimension runs left-right, the next longest runs
up-down, and the camera looks along the shortest, so the front view shows the
broadest face. Pass --view to start from your own rotation matrix instead
(a 3x3 text file, the same format as cellscape's --view).

Usage:
    python cube_views.py structure.cif[.gz] [-o out.png] [--title "Aldolase"]
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import cellscape
from cellscape.cartoon import plot_polygon

# Each face: (name, direction pointing at the camera, direction pointing up on
# screen), both given in the coordinates of the front view. Screen axes are
# x = right, y = up, z = towards the camera, as in cellscape.
FACES = [
    ("front",  (0, 0, 1),  (0, 1, 0)),
    ("left",   (-1, 0, 0), (0, 1, 0)),
    ("right",  (1, 0, 0),  (0, 1, 0)),
    ("back",   (0, 0, -1), (0, 1, 0)),
    ("top",    (0, 1, 0),  (0, 0, -1)),  # looking down, front edge at the bottom
    ("bottom", (0, -1, 0), (0, 0, 1)),   # looking up, front edge at the top
]


def face_matrix(toward_camera, up):
    """Rotation taking front-view coordinates to this face's screen coordinates.

    cellscape rotates atoms as row vectors (``coord @ matrix``), so the matrix
    columns are the new screen x, y and z axes written in the old coordinates.
    """
    z = np.array(toward_camera, dtype=float)
    y = np.array(up, dtype=float)
    x = np.cross(y, z)
    return np.column_stack([x, y, z])


def principal_axes_view(coord):
    """Rotation that puts the longest principal axis on x, then y, then z."""
    centered = coord - coord.mean(axis=0)
    eigenvalues, eigenvectors = np.linalg.eigh(np.cov(centered.T))
    axes = eigenvectors[:, np.argsort(eigenvalues)[::-1]]
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1  # keep it a rotation, not a mirror image
    return axes


def render_face(structure, view_matrix, radius, depth_shading, depth_lines):
    """Outline the structure from one direction and return its styled polygons."""
    structure.set_view_matrix(view_matrix)
    cartoon = structure.outline("residue", radius=radius)
    fig, scratch_axes = plt.subplots()
    cartoon.plot(axes=scratch_axes, do_show=False, color_residues_by="chain",
                 depth_shading=depth_shading, depth_lines=depth_lines)
    plt.close(fig)
    return cartoon._styled_polygons


def polygon_bounds(styled_polygons):
    bounds = np.array([p["polygon"].bounds for p in styled_polygons])
    return bounds[:, 0].min(), bounds[:, 1].min(), bounds[:, 2].max(), bounds[:, 3].max()


def default_output(path):
    name = path.name
    for suffix in (".gz", ".cif", ".mcif", ".mmcif", ".pdb", ".ent"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return path.with_name(f"{name}_cube_views.png")


def make_cube_views(coords, out=None, title=None, view=None, radius=1.5, padding=25,
                    depth_shading=True, depth_lines=True, fig_width=18, dpi=200):
    """Render the six cube-face views of ``coords`` into ``out``. Returns the output path."""
    coords = Path(coords)
    out = Path(out) if out else default_output(coords)
    structure = cellscape.Structure(str(coords), view=False)

    if view is not None:
        base = np.loadtxt(view)
        assert base.shape == (3, 3), "view file must hold a 3x3 matrix"
    else:
        base = principal_axes_view(structure.coord)

    faces = []
    for name, toward_camera, up in FACES:
        matrix = base @ face_matrix(toward_camera, up)
        assert np.isclose(np.linalg.det(matrix), 1.0), "view matrix must be a rotation"
        polygons = render_face(structure, matrix, radius, depth_shading, depth_lines)
        faces.append((name, polygons, polygon_bounds(polygons)))
        print(f"  {name}: {len(polygons)} residues outlined")

    # lay the faces out left to right, centred vertically, at one common scale
    max_height = max(b[3] - b[1] for _, _, b in faces)
    total_width = sum(b[2] - b[0] for _, _, b in faces) + padding * (len(faces) - 1)
    fig_height = fig_width * (max_height * 1.35) / total_width
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    x = 0.0
    for name, polygons, (minx, miny, maxx, maxy) in faces:
        shift = [x - minx, -(miny + maxy) / 2]
        for p in polygons:
            plot_polygon(p["polygon"], axes=ax, translate_pre=shift, facecolor=p["facecolor"],
                         edgecolor=p["edgecolor"], linewidth=p["linewidth"],
                         zorder_mod=p.get("zorder", 0))
        ax.text(x + (maxx - minx) / 2, -max_height / 2 - 0.08 * max_height, name,
                ha="center", va="top", fontsize=14)
        x += (maxx - minx) + padding

    if title:
        ax.text(total_width / 2, max_height / 2 + 0.08 * max_height, title,
                ha="center", va="bottom", fontsize=18)

    ax.set_aspect("equal")
    ax.autoscale()
    ax.axis("off")
    fig.savefig(out, dpi=dpi, bbox_inches="tight", pad_inches=0.1, facecolor="white")
    plt.close(fig)
    print(f"Saved {out}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("coords", help="Structure file (.cif/.pdb, optionally .gz)")
    parser.add_argument("-o", "--out", help="Output PNG (default: <input name>_cube_views.png next to the input)")
    parser.add_argument("--title", help="Title drawn above the views")
    parser.add_argument("--view", help="3x3 rotation matrix file defining the front view (default: principal axes)")
    parser.add_argument("--radius", type=float, default=1.5, help="Atomic radius in angstroms (default: %(default)s)")
    parser.add_argument("--padding", type=float, default=25, help="Gap between views in angstroms (default: %(default)s)")
    parser.add_argument("--no-depth", action="store_true", help="Turn off depth shading and depth line widths")
    parser.add_argument("--dpi", type=int, default=200, help="Output resolution (default: %(default)s)")
    args = parser.parse_args()

    make_cube_views(args.coords, out=args.out, title=args.title, view=args.view,
                    radius=args.radius, padding=args.padding,
                    depth_shading=not args.no_depth, depth_lines=not args.no_depth, dpi=args.dpi)


if __name__ == "__main__":
    main()

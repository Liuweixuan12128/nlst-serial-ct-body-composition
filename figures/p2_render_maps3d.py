"""Paper 2 Figure 7: project the group-level statistical maps onto the 3D paraspinal compartment of one scan.

Every compartment voxel of the example scan (T6-T11) gets the shape-normalised coordinates (u, v) used for the group
maps (pipeline/p2_maps.py) and is assigned the group map value of its cell; the compartment surface is coloured by
probing that volume. Maps: controls' age effect (pp per 10 years) and the case-control level difference (pp).
Writes figures/assets/p2_maps3d_<map>.png (RGBA).

  env -u DISPLAY PYTHONPATH=env/vtk97 VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow \
      env/bin/python figures/p2_render_maps3d.py 100658_T0
"""
import os
import sys

import nibabel as nib
import numpy as np
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
sys.path.insert(0, HERE)
from bodycomp_metrics import CAGE, labels, vertebra_info  # noqa: E402
from render3d import actor, render, surface  # noqa: E402
import vtk  # noqa: E402
from vtk.util import numpy_support as ns  # noqa: E402

NV, NU = 24, 32
SIDES = {"left": ["erector_spinae_left", "transversospinalis_left"],
         "right": ["erector_spinae_right", "transversospinalis_right"]}


def lut(colors, vmin, vmax, n=256):
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list("c", colors)
    t = vtk.vtkLookupTable()
    t.SetNumberOfTableValues(n)
    t.SetRange(vmin, vmax)
    for i in range(n):
        r, g, b, _ = cm(i / (n - 1))
        t.SetTableValue(i, r, g, b, 1.0)
    t.SetNanColor(0.85, 0.85, 0.85, 1.0)
    t.Build()
    return t


def main(cid):
    st = np.load(os.path.join(ROOT, "results", "p2", "map_stats.npz"))
    img = nib.load(os.path.join(ROOT, "data", "nlst_serial", "nifti", f"{cid}.nii.gz"))
    zooms = [float(z) for z in img.header.get_zooms()[:3]]
    seg = os.path.join(ROOT, "data", "seg", "nlst_serial", cid)
    tot = np.asanyarray(nib.load(os.path.join(seg, "total.nii.gz")).dataobj)
    am = np.asanyarray(nib.load(os.path.join(seg, "abdominal_muscles.nii.gz")).dataobj)
    z0, z1 = sorted((int(round(vertebra_info(tot, "T11")[0])), int(round(vertebra_info(tot, "T6")[0]))))
    side_m = {s: labels(am, "abdominal_muscles", n) for s, n in SIDES.items()}
    cord = labels(tot, "total", ["spinal_cord"])
    psm = np.zeros(am.shape, bool)
    for m in side_m.values():
        psm[:, :, z0:z1 + 1] |= m[:, :, z0:z1 + 1]
    out = {"age": st["age_b"], "cc": st["cc_d"]}
    vols = {k: np.full(am.shape, np.nan, np.float32) for k in out}
    for z in range(z0, z1 + 1):
        c = cord[:, :, z]
        xm = np.nonzero(c.any(axis=1))[0].mean() if c.any() else None
        for s, mfull in side_m.items():
            m = mfull[:, :, z]
            if m.sum() < 30:
                continue
            px, py = np.nonzero(m)
            mid = xm if xm is not None else px.mean()
            d = np.abs(px - mid)
            if d.max() - d.min() < 3 or py.max() - py.min() < 3:
                continue
            iu = np.minimum(((d - d.min()) / (d.max() - d.min()) * NU).astype(int), NU - 1)
            iv = np.minimum(((py.max() - py) / (py.max() - py.min()) * NV).astype(int), NV - 1)
            for k, mp in out.items():
                vols[k][px, py, z] = mp[iv, iu]
    # crop for rendering
    zlo, zhi = max(z0 - int(25 / zooms[2]), 0), min(z1 + int(25 / zooms[2]), am.shape[2] - 1)
    sl = (slice(None), slice(None), slice(zlo, zhi + 1))
    bone = labels(tot, "total", CAGE)[sl]
    S_bone = surface(bone, zooms, 0.8)
    S_psm = surface(psm[sl], zooms, 0.9)
    xs, ys = np.nonzero(psm[sl].any(axis=2))
    centre = [(xs.min() + xs.max()) / 2 * zooms[0], (ys.min() + ys.max()) / 2 * zooms[1], (zhi - zlo) / 2 * zooms[2]]
    specs = {"age": (["#f4f3f8", "#cfc9e6", "#9a8fcf", "#6a5bb3", "#4a3aa7", "#2c2170"], 0,
                     float(np.nanpercentile(st["age_b"][st["support"]], 98))),
             "cc": (["#184f95", "#3987e5", "#9ec5f4", "#efeeea", "#f6b896", "#eb6834", "#a8401a"], None, None)}
    v = float(np.nanpercentile(np.abs(st["cc_d"][st["support"]]), 99))
    specs["cc"] = (specs["cc"][0], -v, v)
    for k, (cols, vmin, vmax) in specs.items():
        vol = vols[k][sl]
        # fill a thin shell around the compartment so the surface (slightly outside the mask) samples real values
        filled = vol.copy()
        nanm = np.isnan(filled)
        if nanm.any():
            idx = ndi.distance_transform_edt(nanm, return_distances=False, return_indices=True)
            filled = filled[tuple(idx)]
        im = vtk.vtkImageData()
        im.SetDimensions(*filled.shape)
        im.SetSpacing(*zooms)
        im.GetPointData().SetScalars(ns.numpy_to_vtk(filled.ravel(order="F").astype(np.float32), deep=True,
                                                     array_type=vtk.VTK_FLOAT))
        probe = vtk.vtkProbeFilter()
        probe.SetInputData(S_psm)
        probe.SetSourceData(im)
        probe.Update()
        m = vtk.vtkPolyDataMapper()
        m.SetInputData(probe.GetOutput())
        m.SetLookupTable(lut(cols, vmin, vmax))
        m.SetScalarRange(vmin, vmax)
        m.ScalarVisibilityOn()
        a = vtk.vtkActor()
        a.SetMapper(m)
        a.GetProperty().SetInterpolationToPhong()
        a.GetProperty().SetSpecular(0.1)
        a.GetProperty().SetAmbient(0.25)
        a.GetProperty().SetDiffuse(0.8)
        fn = os.path.join(ROOT, "figures", "assets", f"p2_maps3d_{k}.png")
        render([actor(S_bone, "#f2ecdc", opacity=0.35, spec=0.2), a], centre, -30, 12, (1300, 1700), 1.6, fn)
        print("wrote", fn, vmin, vmax)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "100658_T0")

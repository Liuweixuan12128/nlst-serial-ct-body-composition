"""3D surface renderings of one public CT for Figure 1 (VTK, GPU offscreen via EGL).

Two views of the same patient, identical scale:
  msk  posterior-oblique: skeleton + trunk muscle; muscle inside the T5-T12 span in the tissue red,
       outside the span in a neutral tone (the measured volume is the red block)
  slab the T5-T12 slab alone from anterior-superior: subcutaneous fat, muscle, skeleton and
       visceral/mediastinal fat, from the same slice masks used for the metrics; the top face is the T5 section
Writes figures/assets/render3d_<case>_<view>.png (RGBA, transparent background) and a JSON with the
screen positions of the T4/T8/T12 vertebral-body centroids so Figure 1 can place level ticks exactly.

The project environment keeps VTK 9.3 (needed by TotalSegmentator's fury); VTK 9.7 with EGL lives in env/vtk97:
  env -u DISPLAY PYTHONPATH=env/vtk97 VTK_DEFAULT_OPENGL_WINDOW=vtkEGLRenderWindow \
      env/bin/python figures/render3d.py R01-049
"""
import json
import os
import sys

import numpy as np
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "pipeline"))
sys.path.insert(0, HERE)
from bodycomp_metrics import Case, vertebra_info, labels, CAGE  # noqa: E402
import vtk  # noqa: E402
from vtk.util import numpy_support as ns  # noqa: E402

TISSUE = {"muscle": "#e34948", "sat": "#eda100", "ift": "#2a78d6"}
BONE = "#f2ecdc"
MUSCLE_OUT = "#b9a7a3"
OUT = os.path.join(ROOT, "figures", "assets")


def rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def surface(mask, zooms, sigma_mm=0.9, decimate=0.4):
    """Smoothed iso-surface (0.5) of a binary mask, in mm, RAS axes. Smoothing is at least 0.6 voxel along z
    so thick slices do not show as terraces."""
    f = ndi.gaussian_filter(mask.astype(np.float32),
                            sigma=[sigma_mm / zooms[0], sigma_mm / zooms[1], max(sigma_mm / zooms[2], 0.6)])
    im = vtk.vtkImageData()
    im.SetDimensions(*f.shape)
    im.SetSpacing(*zooms)
    im.GetPointData().SetScalars(ns.numpy_to_vtk(f.ravel(order="F"), deep=True, array_type=vtk.VTK_FLOAT))
    fe = vtk.vtkFlyingEdges3D()
    fe.SetInputData(im)
    fe.SetValue(0, 0.5)
    fe.ComputeNormalsOff()
    fe.Update()
    pd = fe.GetOutput()
    if decimate:
        dc = vtk.vtkQuadricDecimation()
        dc.SetInputData(pd)
        dc.SetTargetReduction(decimate)
        dc.Update()
        pd = dc.GetOutput()
    sm = vtk.vtkWindowedSincPolyDataFilter()
    sm.SetInputData(pd)
    sm.SetNumberOfIterations(20)
    sm.SetPassBand(0.08)
    sm.NormalizeCoordinatesOn()
    sm.Update()
    nr = vtk.vtkPolyDataNormals()
    nr.SetInputData(sm.GetOutput())
    nr.SetFeatureAngle(80)
    nr.SplittingOff()
    nr.ConsistencyOn()
    nr.Update()
    return nr.GetOutput()


def actor(pd, col, opacity=1.0, spec=0.15, amb=0.12, diff=0.85):
    m = vtk.vtkPolyDataMapper()
    m.SetInputData(pd)
    m.ScalarVisibilityOff()
    a = vtk.vtkActor()
    a.SetMapper(m)
    p = a.GetProperty()
    p.SetColor(*rgb(col))
    p.SetOpacity(opacity)
    p.SetSpecular(spec)
    p.SetSpecularPower(25)
    p.SetAmbient(amb)
    p.SetDiffuse(diff)
    p.SetInterpolationToPhong()
    return a


def render(actors, centre, azimuth, elevation, size, zoom, fname, ssao=True):
    ren = vtk.vtkRenderer()
    ren.SetBackground(1, 1, 1)
    ren.SetBackgroundAlpha(0.0)
    for a in actors:
        ren.AddActor(a)
    rw = vtk.vtkRenderWindow()
    rw.SetOffScreenRendering(1)
    rw.SetAlphaBitPlanes(1)
    rw.SetMultiSamples(0)
    rw.AddRenderer(ren)
    rw.SetSize(*size)
    ren.SetUseDepthPeeling(1)
    ren.SetMaximumNumberOfPeels(12)
    ren.SetOcclusionRatio(0.0)
    if ssao:
        ren.UseSSAOOn()
        ren.SetSSAORadius(18.0)
        ren.SetSSAOBias(0.5)
        ren.SetSSAOKernelSize(128)
        ren.SSAOBlurOn()
    ren.UseFXAAOn()
    cam = ren.GetActiveCamera()
    cam.SetFocalPoint(*centre)
    cam.SetPosition(centre[0], centre[1] - 1500, centre[2])   # posterior viewer (+axis1 = anterior)
    cam.SetViewUp(0, 0, 1)
    cam.Azimuth(azimuth)
    cam.Elevation(elevation)
    cam.OrthogonalizeViewUp()
    cam.ParallelProjectionOn()
    ren.ResetCamera()
    cam.Zoom(zoom)
    ren.ResetCameraClippingRange()
    kit = vtk.vtkLightKit()
    kit.SetKeyLightIntensity(0.95)
    kit.SetKeyToFillRatio(2.6)
    kit.SetKeyToHeadRatio(3.0)
    kit.SetKeyToBackRatio(3.5)
    kit.AddLightsToRenderer(ren)
    rw.Render()
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(rw)
    w2i.SetInputBufferTypeToRGBA()
    w2i.ReadFrontBufferOff()
    w2i.Update()
    wr = vtk.vtkPNGWriter()
    wr.SetFileName(fname)
    wr.SetInputConnection(w2i.GetOutputPort())
    wr.Write()
    # screen position (pixels, origin top-left) of given world points
    def to_px(p, _keep=rw):  # the closure keeps the render window alive; without it display coords are 0
        ren.SetWorldPoint(p[0], p[1], p[2], 1.0)
        ren.WorldToDisplay()
        x, y, _ = ren.GetDisplayPoint()
        return float(x), float(size[1] - y)
    return to_px, cam.GetParallelScale()


def main(case_id):
    os.makedirs(OUT, exist_ok=True)
    case = Case(os.path.join(ROOT, "data", "nifti", "public_rg", f"{case_id}.nii.gz"),
                os.path.join(ROOT, "data", "seg", "public_rg", case_id))
    tot = case.segs["total"]
    info = {lv: vertebra_info(tot, lv)[0] for lv in ["T2", "T4", "T5", "T8", "T12", "L1"]}
    z12, z5 = int(round(info["T12"])), int(round(info["T5"]))
    zz = case.zooms[2]
    zlo = max(int(info["T12"] - 40 / zz), 0)
    zhi = min(int(info["T4"] + 55 / zz), tot.shape[2] - 1)
    sl = (slice(None), slice(None), slice(zlo, zhi + 1))
    bone = labels(tot, "total", CAGE + ["scapula_left", "scapula_right", "clavicula_left", "clavicula_right",
                                        "humerus_left", "humerus_right"])[sl]
    mus = case.muscle[sl]
    span = np.zeros(mus.shape[2], bool)
    span[z12 - zlo:z5 - zlo + 1] = True
    mus_in = mus & span[None, None, :]
    mus_out = mus & ~span[None, None, :]
    sat = np.zeros_like(mus)
    ift = np.zeros_like(mus)
    for z in range(z12, z5 + 1):
        m = case.slice_masks(z)
        sat[:, :, z - zlo] = m["sat"]
        ift[:, :, z - zlo] = m["ift"]
    # small islands (noise voxels in the fat window) do not render well
    ift = ndi.binary_opening(ift, structure=np.ones((2, 2, 1), bool))
    lab, n = ndi.label(ift)
    size = ndi.sum(ift, lab, index=np.arange(1, n + 1))
    ift = np.isin(lab, 1 + np.nonzero(size >= 60)[0])
    zooms = case.zooms
    print("building surfaces ...", flush=True)
    S = {k: surface(v, zooms, sigma_mm=s) for k, v, s in [("bone", bone, 0.8), ("mus_in", mus_in, 0.9),
                                                          ("mus_out", mus_out, 0.9), ("sat", sat, 1.6),
                                                          ("ift", ift, 0.9)]}
    body = case.body[sl]
    bx = np.nonzero(body.any(axis=(1, 2)))[0]
    by = np.nonzero(body.any(axis=(0, 2)))[0]
    centre = [(bx.min() + bx.max()) / 2 * zooms[0], (by.min() + by.max()) / 2 * zooms[1],
              (mus.shape[2] - 1) / 2 * zz]
    size = (1400, 1900)
    # vertebral-body centroids (world mm) for level ticks
    cord = labels(tot, "total", ["spinal_cord"])[sl]
    spine_x = np.median(np.nonzero(cord.any(axis=(1, 2)))[0]) * zooms[0]
    spine_y = np.median(np.nonzero(cord.any(axis=(0, 2)))[0]) * zooms[1]
    levels = {lv: [spine_x, spine_y, (info[lv] - zlo) * zz] for lv in ["T4", "T5", "T8", "T12"]}
    meta = {"case": case_id, "size_px": size, "views": {}}
    inspan = lambda m: m & span[None, None, :]
    S["bone_in"] = surface(inspan(bone), zooms, sigma_mm=0.8)
    views = {
        "msk": ([actor(S["bone"], BONE, spec=0.25), actor(S["mus_in"], TISSUE["muscle"]),
                 actor(S["mus_out"], MUSCLE_OUT, opacity=0.45)], -35, 10, (1400, 1900), 1.55),
        "slab": ([actor(S["bone_in"], BONE, spec=0.25), actor(S["mus_in"], TISSUE["muscle"]),
                  actor(S["ift"], TISSUE["ift"], spec=0.2), actor(S["sat"], TISSUE["sat"], spec=0.1)],
                 180 + 30, 38, (1600, 1400), 1.35),
    }
    for name, (acts, az, el, size, zoom) in views.items():
        fn = os.path.join(OUT, f"render3d_{case_id}_{name}.png")
        to_px, scale = render(acts, centre, az, el, size, zoom, fn)
        meta["views"][name] = {"file": os.path.relpath(fn, ROOT), "azimuth": az, "elevation": el, "size_px": size,
                               "mm_per_px": 2 * scale / size[1],
                               "levels_px": {lv: to_px(p) for lv, p in levels.items()}}
        print("wrote", fn, flush=True)
    with open(os.path.join(OUT, f"render3d_{case_id}.json"), "w") as f:
        json.dump(meta, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "R01-049")

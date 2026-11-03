# pvpython --force-offscreen-rendering render_stages.py <dir with mesh_stage_XXX.vtu> <frames dir> [title]
#
# Renders the boundary surface of every adaptive mesh (3D, Fichera domain)
# shaded by log10 of the local mesh size h_T = |T|^{1/3} (light grey =
# coarse, darker grey = fine), with the mesh edges in dark blue. Two PNGs
# per stage: left_XXX.png (whole domain, camera looking into the notch) and
# right_XXX.png (zoom near the re-entrant corner). See README.md for how to
# join them and make an animated GIF.
#
# Created by Ludmil Zikatanov 20261103 with the help of Claude (Anthropic).
import sys, glob, os, math, re
from paraview.simple import *

src, out = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else 'Fichera, P1 Laplace, f = 1'
os.makedirs(out, exist_ok=True)
files = sorted(glob.glob(os.path.join(src, 'mesh_stage_*.vtu')))

def npoints(f):
    with open(f) as fh:
        m = re.search(r'NumberOfPoints="(\d+)" NumberOfCells="(\d+)"', fh.read(400))
    return int(m.group(1)), int(m.group(2))

view = CreateView('RenderView')
view.ViewSize = [1000, 860]
view.OrientationAxesVisibility = 1
view.Background = [1.0, 1.0, 1.0]
try:
    view.UseColorPaletteForBackground = 0
except AttributeError:
    pass
# color range: min over the finest mesh, max over the coarsest
def hrange(f):
    r = XMLUnstructuredGridReader(FileName=[f])
    cs = CellSize(Input=r)
    c = Calculator(Input=cs)
    c.AttributeType = 'Cell Data'
    c.ResultArrayName = 'logh'
    c.Function = 'log10(abs(Volume))/3'
    c.UpdatePipeline()
    a = c.CellData['logh'].GetRange()
    Delete(c); Delete(cs); Delete(r)
    return a[0], a[1]
hmin = hrange(files[-1])[0]
hmax = hrange(files[0])[1]
print('log10 h range', hmin, hmax, flush=True)
lut = GetColorTransferFunction('logh')
lut.RGBPoints = [hmin, 0.80, 0.80, 0.81, hmax, 0.99, 0.99, 0.99]  # grey: fine = dark, coarse = light
lut.ColorSpace = 'RGB'
lut.RescaleTransferFunction(hmin, hmax)
lut.AutomaticRescaleRangeMode = 'Never'

text = Text()
tdisp = Show(text, view)
tdisp.FontSize = 22
tdisp.Color = [0, 0, 0]
tdisp.WindowLocation = 'Upper Left Corner'

prev = []
for k, f in enumerate(files):
    for p in prev:
        Delete(p)
    nv, ns = npoints(f)
    r = XMLUnstructuredGridReader(FileName=[f])
    cs = CellSize(Input=r)
    cs.ComputeVertexCount = 0; cs.ComputeLength = 0; cs.ComputeArea = 0
    calc = Calculator(Input=cs)
    calc.AttributeType = 'Cell Data'
    calc.ResultArrayName = 'logh'
    calc.Function = 'log10(abs(Volume))/3'
    surf = ExtractSurface(Input=calc)
    d = Show(surf, view)
    d.SetRepresentationType('Surface With Edges')
    ColorBy(d, ('CELLS', 'logh'))
    d.LookupTable = lut
    lut.AutomaticRescaleRangeMode = 'Never'
    lut.RescaleTransferFunction(hmin, hmax)
    d.EdgeColor = [0.0, 0.0, 0.45]
    d.Ambient = 0.6   # lighting: keep faces turned away from the light light grey
    d.Diffuse = 0.45
    d.LineWidth = 0.6
    if k == 0:
        d.SetScalarBarVisibility(view, True)
        bar = GetScalarBar(lut, view)
        bar.Title = 'log10 h_T'
        bar.ComponentTitle = ''
        bar.TitleColor = [0, 0, 0]
        bar.LabelColor = [0, 0, 0]
        bar.TitleFontSize = 18
        bar.LabelFontSize = 16
    text.Text = title + '\nstage %d of SOLVE-ESTIMATE-MARK-REFINE\n%d vertices, %d tetrahedra' % (k, nv, ns)
    # left: whole domain, with edges
    d.SetRepresentationType('Surface With Edges')
    tdisp.Visibility = 1
    text.Text = title + '\nstage %d of SOLVE-ESTIMATE-MARK-REFINE\n%d vertices, %d tetrahedra' % (k, nv, ns)
    view.CameraPosition = [4.2, -3.6, 3.4]
    view.CameraFocalPoint = [0.0, 0.0, 0.0]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.CameraParallelScale = 1.9
    view.CameraViewAngle = 30
    Render(view)
    SaveScreenshot(os.path.join(out, 'left_%03d.png' % k), view, ImageResolution=[1000, 860])
    # right: zoom on the re-entrant corner, with edges
    d.SetRepresentationType('Surface With Edges')
    text.Text = 'zoom near the re-entrant corner (0,0,0)'
    view.CameraPosition = [1.45, -1.25, 1.2]
    view.CameraFocalPoint = [0.12, -0.12, 0.12]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    Render(view)
    SaveScreenshot(os.path.join(out, 'right_%03d.png' % k), view, ImageResolution=[1000, 860])
    prev = [surf, calc, cs, r]
    print('frame', k, nv, ns, flush=True)

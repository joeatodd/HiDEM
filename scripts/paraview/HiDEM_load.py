try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

HiDEM_data = GetActiveSource()

renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

glyph1 = Glyph(Input=HiDEM_data,
    GlyphType='2D Glyph')
glyph1.GlyphMode = 'All Points'
glyph1.GlyphType.GlyphType = 'Vertex'

glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 645.2785175323487
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 32.263925876617435
glyph1Display.SetScaleArray = [None, '']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = [None, '']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

Hide(HiDEM_data, renderView1)

# hide data in view
Hide(glyph1, renderView1)

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1)

# Properties modified on glyph1Display
glyph1Display.RenderPointsAsSpheres = 1

# Properties modified on glyph1Display
glyph1Display.PointSize = 8.0

renderView1.Update()

# hide data in view
#Hide(jYR00, spreadSheetView1)

SetActiveView(renderView1)

# reset view to fit data
renderView1.ResetCamera()
renderView1.Update()

Render()

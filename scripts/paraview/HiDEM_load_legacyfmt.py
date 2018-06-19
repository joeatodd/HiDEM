try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

HiDEM_data = GetActiveSource()

renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

# Properties modified on HiDEM_data
HiDEM_data.HaveHeaders = 0
HiDEM_data.FieldDelimiterCharacters = ' '
HiDEM_data.MergeConsecutiveDelimiters = 1

renderView1.Update()
# create a new 'Table To Points'
tableToPoints1 = TableToPoints(Input=HiDEM_data)
tableToPoints1.XColumn = 'Field 0'
tableToPoints1.YColumn = 'Field 1'
tableToPoints1.ZColumn = 'Field 2'

# show data in view
#tableToPoints1Display = Show(tableToPoints1, spreadSheetView1)
tableToPoints1Display = Show(tableToPoints1, renderView1)

# hide data in view
#Hide(jYR00, spreadSheetView1)

SetActiveView(renderView1)

tableToPoints1Display.Representation = 'Surface'
tableToPoints1Display.ColorArrayName = [None, '']
tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
tableToPoints1Display.SelectOrientationVectors = 'None'
tableToPoints1Display.ScaleFactor = 647.6596292
tableToPoints1Display.SelectScaleArray = 'None'
tableToPoints1Display.GlyphType = 'Arrow'
tableToPoints1Display.GlyphTableIndexArray = 'None'
tableToPoints1Display.DataAxesGrid = 'GridAxesRepresentation'
tableToPoints1Display.PolarAxes = 'PolarAxesRepresentation'
tableToPoints1Display.GaussianRadius = 323.8298146
tableToPoints1Display.SetScaleArray = [None, '']
tableToPoints1Display.ScaleTransferFunction = 'PiecewiseFunction'
tableToPoints1Display.OpacityArray = [None, '']
tableToPoints1Display.OpacityTransferFunction = 'PiecewiseFunction'


# Properties modified on renderView1
#renderView1.EnableOSPRay = 1

# Properties modified on tableToPoints1Display
tableToPoints1Display.PointSize = 8.0
tableToPoints1Display.RenderPointsAsSpheres = 1

# reset view to fit data
renderView1.ResetCamera()
renderView1.Update()

Render()

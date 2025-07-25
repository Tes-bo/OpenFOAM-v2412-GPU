#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

useGroup = True
pieces = []

# Load all pieces
for worldName in ['top', 'slab1', 'slab2', 'slab3', 'slab4']:
    # Could also touch world/world.foam into existence
    loadedDir = './' + worldName + '/' + worldName + '.foam'
    piece = OpenFOAMReader(registrationName=worldName, FileName=loadedDir)
    piece.SkipZeroTime = 0
    piece.Decomposepolyhedra = 0
    piece.MeshRegions = ['internalMesh']
    piece.CellArrays = ['T']
    pieces.append(piece)

# Set display for all pieces
for piece in pieces:
    display = Show(piece, renderView1, 'UnstructuredGridRepresentation')
    display.SetRepresentationType('Surface With Edges')
    ColorBy(display, ('CELLS', 'T'))
    # Show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)

# ----
# Group pieces
if useGroup:
    # Hide data in view
    group1 = GroupDatasets(registrationName='GroupDatasets1', Input=pieces)

    # show data in view
    display = Show(group1, renderView1, 'UnstructuredGridRepresentation')

    display.SetRepresentationType('Surface With Edges')
    ColorBy(display, ('CELLS', 'T'))
    # Show color bar/color legend
    display.SetScalarBarVisibility(renderView1, True)

    for piece in pieces:
        Hide(piece, renderView1)

# ----
# get color transfer function/color map for 'T'
tLUT = GetColorTransferFunction('T')

# Rescale transfer function
tLUT.RescaleTransferFunction(0.0, 1.2)

# get opacity transfer function/opacity map for 'T'
tPWF = GetOpacityTransferFunction('T')

# Rescale transfer function
tPWF.RescaleTransferFunction(0.0, 1.2)

# reset view to fit data
renderView1.ResetCamera(False)

# update the view to ensure updated data information
renderView1.Update()

#================================================================

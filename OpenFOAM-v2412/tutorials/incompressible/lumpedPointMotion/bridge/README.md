Setting up the FSI linear controllers can be easy or difficult,
depending on the complexity of the structure and the source of the
input data.

If you already have an FEA model with nodes and beam elements, it will
be fairly straightforward to instrument your structure. If you just
have the control points but need to define the connectivity yourself,
it will obviously be more work.

In this case it can help to import your raw points as a csv table in
paraview and use the TableToPoints filter to visualize their locations
and overlay this with the surface geometry. With this you can piece
together the connectivity, which specifies the motion _controllors_.

To debug the setup, within `steady/`

Copy the geometry
```
./Allrun.init
```
Verify connectivity
```
lumpedPointZones -dry-run
```

This will generate a `state.vtp` file with the lumped points, connected
as per the controller description.


Next generate the mesh (eg, with snappyHexMesh). For example,
```
./Allrun.pre
```

Test the mapping
```
lumpedPointZones   # serial or parallel

paraview lumpedPointZones.vtp
```
Inspect the nearest/next nearest and weighting.
Adjust the controllers definitions if required.


For setup, it is often helpful if you have some predefined structural
response data that can be used for testing.

Check the quality of response data by visualizing how it affects the
movement of the points:
```
lumpedPointMovement {options} -dry-run response.txt

paraview state.vtk.series
```

You can add visualization options such as `-scale` or
`-visual-length`. If the there are many time points in the response
data, use `-span` to skip over some of them.


The `lumpedPointMovement` command can be used as above, but without
the `-dry-run` option. This will extract the patch surface associated
with the controllers and generate corresponding surface files in VTK
format.
```
lumpedPointMovement {options} response.txt

paraview geom.vtp.series
```

Using a larger scale factor (eg, `-scale 10`) can help highlight
potential interpolation problems.


## Additional Notes

The `files/polynomials.dict` represent a vague approximation of
measurement data, but are only intended as a compacter means of
representing the movement. They can be used with the accompanying
`code/polynomial-motion.C` to act as a response slave, or to
pre-generate a response file. For example,

```
code/polynomial-motion -output response.txt -deltaT 0.001 -nTimes 5001 files/polynomials.dict
```

To query the values at an individual time:
```
code/polynomial-motion -query -time 0.5 files/polynomials.dict
```

## Preliminary layout for OpenFOAM data within ADIOS

* This is a work-in-progress (proof-of-concept) for reading/writing
  OpenFOAM data with ADIOS.
* For development purposes, a function object is used for the implementation framework.
* Restart is with a 'fast-forward' principle: OpenFOAM is used for the Time=0
  construction of meshes and fields, fields are overwritten with their
  values read from the adios file.
* Supports multiple regions, multiple clouds per region.
* Restarting with mesh changing is not yet fully implemented.

---

2019-12-20

---

* Storage directory:  "adiosData/"

* Each "bp" (OpenFOAM binary packed) adios file is restricted to a
  ***single*** time-step/iteration.
  This makes for simple and efficient handling.

* In rare cases are values stored in global arrays with offsets.
  - Single values per-processor (eg, time index).
  - Cloud parcel data (currently).

  In all other cases, data are stored with their local dimensions only.

* Data stored as "unsigned byte" generally represent binary content
  that has been serialized via the OpenFOAM `Ostream` and are targeted
  for use by an `Istream` when reading. This binary content should be
  largely identical to a normal OpenFOAM field-file, but without the
  file-header.

---

### General Attributes

These attributes provide assistance when reading the data files.
All entries are considered mandatory.


| type    | name                        | example
|---------|-----------------------------|--------
| _any_   | /constant/...               | _reserved_
| _any_   | /system/...                 | _reserved_
| int     | /openfoam/api               | 1906
| string  | /openfoam/patch             | "190401"
| string  | /openfoam/endian            | "LSB"
| int     | /openfoam/label             | 32
| int     | /openfoam/scalar            | 64
| int     | /openfoam/nProcs            | 4
| int     | /openfoam/nRegions          | 2
| string[]| /openfoam/regions           | {"region0", "solid"}


The number of regions can either be obtained directly from the
`/openfoam/nRegions` attribute, or aternatively from the list length
of the `/openfoam/regions` attribute. No particular sort order is
specified for the region names.


### General Variables

Coordinated data for all mesh, fields and cloud information.

| type    | name                        | comment
|---------|-----------------------------|-------------
| _any_   | /constant/...               | _reserved_
| _any_   | /system/...                 | _reserved_


### Time Attributes

Since each file is restricted to a ***single*** time-step/iteration,
which inherently identical across all processes, the time management
values are tracked as attributes rather than as variables. This
additionally simplifies later post-processing, since the attribute
values are available directly from the meta-data.


| type    | name            | comment
|---------|-----------------|-------------
| int     | /time/index     | iteration value
| double  | /time/value     | time-value
| double  | /time/deltaT    | current time-step value
| double  | /time/deltaT0   | previous time-step value
| double  | /time/faces     | time-value of mesh topology (per region)
| double  | /time/points    | time-value of mesh points (per region)

The order of the `/time/faces` and `/time/points` attributes must
correspond exactly to the ordering used by the `/openfoam/regions` attribute.

The attribute `/time/value` can be considered to be a global time value
for all fields and clouds contained within the file.
When the `/time/faces` value for a particular region matches that of
`/time/value`, the file will also contain a mesh description for that
region.
If only the `/time/points` value for a particular region matches, this
indicates mesh motion (but with the same topology) and the file will
correspondingly contain the updated points for that region.

This mechanism can be used to locate an appropriate mesh description from
previously saved files.
The following pseudocode illustrates part of the logic:

    if /time/value == /time/faces/REGION
        # current file contains updated topology (points, faces, etc)
        REGION/polyMesh/points
        REGION/polyMesh/faces/...
        ...

    elif /time/value == /time/points/REGION
        # current file contains updated mesh points
        REGION/polyMesh/points

    else
        # mesh points, faces in an older file
        # - compare to currently loaded meshes to see if update is needed


Since the adios files are generally saved with names corresponding to
their time state, it is possible to identify the file containing a
particular mesh/time-state without searching the file contents.


### Region Attributes

Basic information that applies to mesh and field information.

| type    | name                        | comment
|---------|-----------------------------|--------
| int     | \<regionName\>/nPatches     | number of patches
| string[]| \<regionName\>/patch-names  | patch names
| string[]| \<regionName\>/patch-types  | patch types

The number of patches used in the mesh may be larger than the number
of patch names and types attributes since these attributes do not
include any processor patches.


#### Example

| type    | name                        | example
|---------|-----------------------------|--------
| int     | region0/nPatches            | 3
| string[]| region0/patch-names         | {"movingWall", "fixedWalls", "frontAndBack"}
| string[]| region0/patch-types         | {"wall", "wall", "empty" }
|         | solid/...                   |


### Meshes

These variables are only available when the *updateMesh* attribute is
also true.

#### Global Array Variables


| type    | name                                | comment
|---------|-------------------------------------|-------------
| int     | \<regionName\>/polyMesh/nPoints     | per processor
| int     | \<regionName\>/polyMesh/nCells      | per processor
| int     | \<regionName\>/polyMesh/nFaces      | per processor
| int     | \<regionName\>/polyMesh/nInternalFaces | per processor


#### Local Variables

| type    | name                                  | comment
|---------|---------------------------------------|-------------
| double  | \<regionName\>/polyMesh/points        | {nPoints x 3}
| int     | \<regionName\>/polyMesh/faces/indices | indices for compact faceList
| int     | \<regionName\>/polyMesh/faces/content | content for compact faceList
| int     | \<regionName\>/polyMesh/owner         | face owners
| int     | \<regionName\>/polyMesh/neighbour     | face neighbours
| byte    | \<regionName\>/polyMesh/boundary      | boundary mesh information


### Fields


#### Local Variables

| type    | name                                 | comment
|---------|--------------------------------------|-------------
| double  | \<regionName\>/field/U               | primitive field content
| byte    | \<regionName\>/field/U/boundaryField | internalField - primitive field content
| double  | \<regionName\>/field/p               | internalField - primitive field content
| byte    | \<regionName\>/field/p/boundaryField | as per OPENFOAM binary file content
|         | ...


#### Attributes

| type    | name                                | example
|---------|-------------------------------------|-------------
| string  | \<regionName\>/field/U/class        | "volVectorField"
| string[]| \<regionName\>/field/U/patch-types  | {"fixedValue", "fixedValue", "empty"}
| string  | \<regionName\>/field/p/class        | "volScalarField"
| string[]| \<regionName\>/field/p/patch-types  | { "zeroGradient", "zeroGradient", "empty" }
|         | ...


### Clouds

These entries are only available when the region also has any active
cloud information. A missing value of `nClouds` is equivalent to `nClouds=0`.

The cloud contents are written as their IOFields.

#### Region Cloud Attributes

| type    | name                        | comment
|---------|-----------------------------|-------------
| int     | \<regionName\>/nClouds      | number of associated clouds in region
| string[]| \<regionName\>/clouds       | {"name0", "name1", ...}


#### Cloud Variables (Globally-addressable)

| type    | name                           | comment
|---------|--------------------------------|-------------
| double  | \<regionName\>/cloud/\<cloudName\>/position | cloud x,y,z positions


#### Cloud Attributes

| type    | name                           | example
|---------|--------------------------------|-------------
| string  | \<regionName\>/cloud/\<cloudName\>/class    | "Cloud\<basicKinematicCollidingParcel\>"
| int     | \<regionName\>/cloud/\<cloudName\>/nParcels | total count (sum of corresponding variable)


---

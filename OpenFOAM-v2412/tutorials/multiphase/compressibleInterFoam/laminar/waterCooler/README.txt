Simple testcase to show multi-world explicit coupling. It uses 'solidFoam'
(the solid functionality from chtMultiRegionFoam) and 'compressibleInterFoam'
to demonstrate cooling a solid block with water. WIP.


World solid
-----------
This is running in a single directory 'solid'. Its local patch 'top'
communicates with the 'fluid' world, patch 'bottom'. It runs solidFoam
to simulate a heated solid block. Note that the physical settings probably
do not make much sense.


World fluid
-----------
This is running in the 'fluid' directory. The local patch 'bottom'
communicates with the 'top' patch of the solid world.

Changes to the original tutorial ('depthCharge2D'):
- changed to a fixed timestep.
- added the 'syncObjects' functionObject. This transfers the data from
one world to the other, once every timestep.
- changed the system/blockMeshDict to define the coupling on the patch:
    bottom
    {
        type mappedWall;

        // What to sample:
        sampleMode      nearestPatchFace;

        // Simulation world to sample
        sampleWorld     solid;

        // Region to sample
        sampleRegion    region0;

        // If sampleMode is nearestPatchFace : patch to find faces of
        samplePatch     top;

        // Use database to get data from (one-way or loose coupling in
        // combination with functionObject)
        sampleDatabase    true;

        faces
        (
            (1 5 4 0)
        );
    }


A similar blockMeshDict boundary definition was used on the solid side:
    top
    {
        type mappedWall;

        // What to sample:
        sampleMode      nearestPatchFace;

        // Simulation world to sample
        sampleWorld     fluid;

        // Region to sample
        sampleRegion    region0;

        // If sampleMode is nearestPatchFace : patch to find faces of
        samplePatch     bottom;

        // Use database to get data from (one-way or loose coupling in
        // combination with functionObject)
        sampleDatabase    true;
    }


Note that if the 'sampleDatabase' is not mentioned the coupling will be done
from boundary condition to boundary condition. This is a closer coupling but
does mean that both solvers have to run the same number of correctors.

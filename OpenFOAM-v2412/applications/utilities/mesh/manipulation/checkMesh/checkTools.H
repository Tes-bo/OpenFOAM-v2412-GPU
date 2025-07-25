#include "scalar.H"
#include "indirectPrimitivePatch.H"

namespace Foam
{
    // Forward Declarations
    class polyMesh;
    class pointSet;
    class faceSet;
    class cellSet;
    class fileName;
    class polyMesh;
    class coordSetWriter;
    class surfaceWriter;

    void printMeshStats(const polyMesh& mesh, const bool allTopology);

    //- Generate merged surface on master and write. Needs input patch
    //  to be of mesh faces.
    void mergeAndWrite
    (
        const polyMesh& mesh,
        surfaceWriter& writer,
        const word& name,
        const indirectPrimitivePatch& setPatch,
        const fileName& outputDir
    );

    //- Write vtk representation of (assembled) faceSet to surface file in
    //  postProcessing/ directory
    void mergeAndWrite(surfaceWriter& writer, const faceSet& set);

    //- Write vtk representation of (assembled) cellSet to surface file in
    //  postProcessing/ directory
    void mergeAndWrite(surfaceWriter& writer, const cellSet& set);

    //- Write vtk representation of (assembled) pointSet to 'set' file in
    //  postProcessing/ directory
    void mergeAndWrite(coordSetWriter& writer, const pointSet& set);
}


// ************************************************************************* //

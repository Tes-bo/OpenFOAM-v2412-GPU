#include "polyMesh.H"

namespace Foam
{
    class surfaceWriter;

    label checkMeshQuality
    (
        const polyMesh& mesh,
        const dictionary& dict,
        autoPtr<surfaceWriter>& writer
    );
}

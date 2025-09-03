#include "label.H"
#include "HashSet.H"
#include "labelVector.H"

namespace Foam
{
    // Forward Declarations
    class polyMesh;
    class wedgePolyPatch;
    class coordSetWriter;
    class surfaceWriter;

    label findOppositeWedge(const polyMesh&, const wedgePolyPatch&);

    //- Check wedge orientation
    bool checkWedges
    (
        const polyMesh&,
        const bool report,
        const Vector<label>&,
        labelHashSet*
    );

    //- Check 0th vertex on coupled faces
    bool checkCoupledPoints(const polyMesh&, const bool report, labelHashSet*);

    //- Collect AMI weights to master and write
    void collectAndWriteAMIWeights
    (
        const polyMesh& mesh,
        surfaceWriter& wr,
        const fileName& fName,
        const scalarField& weights,
        const faceList& localFaces,
        const labelList& meshPoints,
        const Map<label>& meshPointMap,

        // Collect geometry
        faceList& mergedFaces,
        pointField& mergedPoints,
        autoPtr<globalIndex>& globalFaces,
        autoPtr<globalIndex>& globalPoints
    );

    label checkGeometry
    (
        const polyMesh& mesh,
        const bool allGeometry,
        autoPtr<surfaceWriter>& surfWriter,
        autoPtr<coordSetWriter>& setWriter
    );
}

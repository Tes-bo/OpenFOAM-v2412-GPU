#include "createMesh.H"

dictionary filmDict;

IOobject io
(
    "surfaceFilmProperties",
    mesh.time().constant(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::NO_WRITE,
    IOobject::NO_REGISTER
);

if (io.typeHeaderOk<IOdictionary>())
{
    IOdictionary propDict(io);

    filmDict = std::move(propDict);

    const word filmRegionName = filmDict.get<word>("region");

    fvMesh filmMesh
    (
        IOobject
        (
            filmRegionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
}

// ************************************************************************* //

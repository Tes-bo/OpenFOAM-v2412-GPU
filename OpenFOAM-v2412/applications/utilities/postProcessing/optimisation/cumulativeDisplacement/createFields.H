    pointField points0
    (
        pointIOField
        (
            IOobject
            (
                "points",
                mesh.time().constant(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        )
    );

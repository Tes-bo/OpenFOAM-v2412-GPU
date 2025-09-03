{
    IOobject io
    (
        "points",
        runTime.timeName(),
        mesh.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (io.typeHeaderOk<pointIOField>(true, false))
    {
        Info<< "new points available" << endl;
        // Read new points
        io.readOpt(IOobject::MUST_READ);
        mesh.movePoints(pointIOField(io));
    }
}

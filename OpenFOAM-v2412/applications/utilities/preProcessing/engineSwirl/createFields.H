Info<< "Reading combustion properties\n" << endl;

IOdictionary engineGeometry
(
    IOobject
    (
        "engineGeometry",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

vector swirlAxis(engineGeometry.get<vector>("swirlAxis"));

vector swirlCenter(engineGeometry.get<vector>("swirlCenter"));

dimensionedScalar swirlRPMRatio("swirlRPMRatio", engineGeometry);

dimensionedScalar swirlProfile("swirlProfile", engineGeometry);

dimensionedScalar bore("bore", dimLength, engineGeometry);

dimensionedScalar rpm("rpm", dimless/dimTime, engineGeometry);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

vector zT = swirlAxis;
vector yT = vector(0, zT.z(), -zT.y());
vector xT = vector
(
    zT.y()*zT.y() + zT.z()*zT.z(),
    -zT.x()*zT.y(),
    -zT.x()*zT.z()
);

// if swirl is around (1, 0, 0) we have to find another transformation
if (mag(yT) < SMALL)
{
    yT = vector(zT.y(), -zT.x(), 0);
    xT = vector(-zT.x()*zT.z(), -zT.y()*zT.z(), zT.x()*zT.x() + zT.y()*zT.y());
}

//swirlAxis doesn't have to be of unit length.
xT.normalise();
yT.normalise();
zT.normalise();

IOdictionary dict
(
    IOobject
    (
        "createBoxTurbDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ
    )
);

// Extents in x, y, z directions
const vector L(dict.get<vector>("L"));

// Number of cells in x, y, z directions
const Vector<label> N(dict.get<Vector<label>>("N"));

// Wave number vs energy profile
autoPtr<Function1<scalar>> Ek(Function1<scalar>::New("Ek", dict, &runTime));

// Number of modes
const label nModes = dict.get<label>("nModes");

// Mesh spacing in x, y and z directions
const vector delta
(
    L.x()/scalar(N.x()),
    L.y()/scalar(N.y()),
    L.z()/scalar(N.z())
);

Random rndGen(1234);

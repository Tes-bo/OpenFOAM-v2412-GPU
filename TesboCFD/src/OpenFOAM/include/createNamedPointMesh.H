Foam::Info
    << "Create pointMesh for time = "
    << runTime.timeName() << Foam::nl << Foam::endl;

// Register pointMesh on the database
const Foam::pointMesh& pMesh = pointMesh::New
(
    mesh,
    Foam::IOobject::READ_IF_PRESENT
);

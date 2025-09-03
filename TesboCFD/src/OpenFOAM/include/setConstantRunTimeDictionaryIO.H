IOobject dictIO = IOobject::selectIO
(
    IOobject
    (
        dictName,
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    ),
    args.getOrDefault<fileName>("dict", "")
);

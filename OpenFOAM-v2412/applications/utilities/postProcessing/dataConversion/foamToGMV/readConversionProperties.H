IOdictionary conversionProperties
(
    IOobject
    (
        "conversionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


scalar startTime
(
    conversionProperties.get<scalar>("startTime")
);

word vComp
(
    conversionProperties.get<word>("vector")
);

word format
(
    conversionProperties.get<word>("format")
);

if ((format != "ascii") && (format != "ieeei4r8"))
{
    FatalError
        << "format type: " << format << " unknown."
        << " Valid options are: ascii ieeei4r8"
        << abort(FatalError);

}

word cells
(
    conversionProperties.get<word>("cells")
);

if
(
    (cells != "hex")
)
{
    FatalError
        << "cells type: " << cells << " unknown."
        << abort(FatalError);
}

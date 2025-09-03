/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Read time index from ../uniform/time, but treat 0 and constant specially
    or simply increment from the '-index' option if it was supplied

Note
    Does a 'continue' out of the loop if not possible to obtain an index

\*---------------------------------------------------------------------------*/

label timeIndex = 0;
{
    bool goodTimeIndex = true;

    if (doConsecutive)
    {
        timeIndex = indexingNumber++;
    }
    else if
    (
        runTime.timeName() != runTime.constant()
     && runTime.timeName() != "0"
    )
    {
        IOobject io
        (
            "time",
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        if (io.typeHeaderOk<IOdictionary>(true, false))
        {
            io.readOpt(IOobject::MUST_READ_IF_MODIFIED);
            IOdictionary timeObject(io);

            timeObject.readEntry("index", timeIndex);
        }
        else
        {
            goodTimeIndex = false;
            Info<< "skip ... missing file: " << io.objectRelPath() << endl;
        }
    }

    if (!goodTimeIndex)
    {
        continue;
    }
}


// ************************************************************************* //

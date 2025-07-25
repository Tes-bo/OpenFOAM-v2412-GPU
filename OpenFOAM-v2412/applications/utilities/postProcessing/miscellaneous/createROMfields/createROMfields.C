/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    createROMfields

Group
    grpPostProcessingUtilities

Description
    Create fields using reduced-order modelling (ROM) data
    at specific time instants without requiring any CFD computations.

Usage
    Minimal example by using \c system/ROMfieldsDict:
    \verbatim
        // Mandatory entries
        ROMmodel        <word>;

        // Inherited entries
        // See DMD.H for the 'DMD' ROMmodel
        ...
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type | Reqd | Deflt
      ROMmodel  | Type of reduced-order model    | word | yes  | -
    \endtable

    Options for the \c ROMmodel entry:
    \verbatim
      DMD    | Streaming total dynamic mode decomposition
    \endverbatim

    The inherited entries are elaborated in:
      - \link ROMmodel.H \endlink
      - \link DMD.H \endlink

Note
  - The quality of results depends on the capabilities of the underlying
    reduced-order model, and the quality of the input data.
  - Warning: Reduced-order modelling is an active research area at the time of
    writing; therefore, there could be cases whereat oddities can be seen.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "fvCFD.H"
#include "ROMmodels/ROMmodel/ROMmodel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create fields using reduced-order modelling (ROM) data at specific "
        "time instants without requiring any CFD computations."
    );

    argList::addOption
    (
        "dict",
        "file",
        "Alternative dictionary for ROMfieldsDict"
    );

    // No -constant, no special treatment for 0/
    timeSelector::addOptions(false);

    // Remove treatments unnecessary for this utility
    argList::noFunctionObjects();
    argList::removeOption("noZero");
    argList::removeOption("world");


    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("ROMfieldsDict");
    #include "setSystemRunTimeDictionaryIO.H"
    Info<< "Reading " << dictIO.name() << nl << endl;
    IOdictionary dict(dictIO);

    instantList times = timeSelector::select0(runTime, args);
    if (times.empty())
    {
        FatalErrorInFunction
            << "No times selected." << nl
            << exit(FatalError);
    }

    #include "createNamedMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    auto ROMptr = ROMmodel::New(runTime, mesh, dict, times);

    ROMptr->read(dict);

    ROMptr->createAndWrite();


    Info<< nl;
    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    chemkinToFoam

Group
    grpSurfaceUtilities

Description
    Convert CHEMKINIII thermodynamics and reaction data files into
    OpenFOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "chemkinReader.H"
#include "OFstream.H"
#include "StringStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // More precision (for output of JANAF coefficients)
    IOstream::minPrecision(10);

    argList::addNote
    (
        "Convert CHEMKINIII thermodynamics and reaction data files into"
        " OpenFOAM format."
    );

    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects

    argList::addArgument("CHEMKINFile");
    argList::addArgument("CHEMKINThermodynamicsFile");
    argList::addArgument("CHEMKINTransport");
    argList::addArgument("FOAMChemistryFile");
    argList::addArgument("FOAMThermodynamicsFile");

    argList::addBoolOption
    (
        "newFormat",
        "Read Chemkin thermo file in new format"
    );

    argList args(argc, argv);

    const bool newFormat = args.found("newFormat");

    speciesTable species;

    chemkinReader cr
    (
        species,
        args.get<fileName>(1),  // chemkin fileName
        args.get<fileName>(3),  // thermo fileName
        args.get<fileName>(2),  // transport fileName
        newFormat
    );

    {
        // output: reactions file
        OFstream reactionsFile(args.get<fileName>(4));

        reactionsFile.writeEntry("elements", cr.elementNames()) << nl;
        reactionsFile.writeEntry("species", cr.species()) << nl;

        cr.reactions().write(reactionsFile);
    }

    // Temporary hack to splice the specie composition data into the thermo file
    // pending complete integration into the thermodynamics structure

    OCharStream os;
    cr.speciesThermo().write(os);

    ISpanStream is(os.view());
    dictionary thermoDict(is);

    // Add elements
    for (entry& dEntry : thermoDict)
    {
        const word& speciesName = dEntry.keyword();
        dictionary& speciesDict = dEntry.dict();

        dictionary elemDict("elements");

        for (const specieElement& elem : cr.specieComposition()[speciesName])
        {
            elemDict.add(elem.name(), elem.nAtoms());
        }

        speciesDict.add("elements", elemDict);
    }

    // output: thermo file

    thermoDict.write(OFstream(args.get<fileName>(5))(), false);


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

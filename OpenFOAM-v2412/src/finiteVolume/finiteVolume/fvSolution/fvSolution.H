/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

Class
    Foam::fvSolution

Description
    Selector class for finite volume solution solution.
    fvMesh is derived from fvSolution so that all fields have access to the
    fvSolution from the mesh reference they hold.

    The optional fallback dictionary content for constructors is used
    when a file is missing or for a NO_READ, with a null pointer being
    treated like an empty dictionary.

SeeAlso
    Foam::solution

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvSolution_H
#define Foam_fvSolution_H

#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvSolution Declaration
\*---------------------------------------------------------------------------*/

class fvSolution
:
    public solution
{
public:

    //- No copy construct
    fvSolution(const fvSolution&) = delete;

    //- No copy assignment
    void operator=(const fvSolution&) = delete;


    // Constructors

        //- Construct for objectRegistry, readOption, (system) dictionary name.
        fvSolution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, rOpt, dictName, fallback)
        {}

        //- Construct for objectRegistry, dictionary name and optional
        //- fallback dictionary content (for a NO_READ or missing file)
        //  A null dictionary pointer is treated like an empty dictionary.
        fvSolution
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, dictName, fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("fvSolution").
        fvSolution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, rOpt, "fvSolution", fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("fvSolution") and
        //- fallback dictionary content.
        fvSolution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary& dict
        )
        :
            solution(obr, rOpt, "fvSolution", &dict)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("fvSolution").
        //  Uses the readOption from the registry.
        explicit fvSolution
        (
            const objectRegistry& obr,
            const dictionary* fallback = nullptr
        )
        :
            fvSolution(obr, "fvSolution", fallback)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("fvSolution") and optional contents.
        //  Uses the readOption from the registry.
        fvSolution(const objectRegistry& obr, const dictionary& dict)
        :
            fvSolution(obr, "fvSolution", &dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

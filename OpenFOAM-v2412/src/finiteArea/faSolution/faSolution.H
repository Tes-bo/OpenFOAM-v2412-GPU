/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::faSolution

Description
    Selector class for finite area solution.
    faMesh is derived from faSolution so that all fields have access to the
    faSolution from the mesh reference they hold.

    The optional fallback dictionary content for constructors is used
    when a file is missing or for a NO_READ, with a null pointer being
    treated like an empty dictionary.

SeeAlso
    Foam::solution

\*---------------------------------------------------------------------------*/

#ifndef Foam_faSolution_H
#define Foam_faSolution_H

#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faSolution Declaration
\*---------------------------------------------------------------------------*/

class faSolution
:
    public solution
{
public:

    // Generated Methods

        //- No copy construct
        faSolution(const faSolution&) = delete;

        //- No copy assignment
        void operator=(const faSolution&) = delete;


    // Constructors

        //- Construct for objectRegistry, readOption, (system) dictionary name.
        faSolution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, rOpt, dictName, fallback)
        {}

        //- Construct for objectRegistry, (system) dictionary name
        //- using the readOption from the registry.
        faSolution
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, dictName, fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("faSolution").
        faSolution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, rOpt, "faSolution", fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("faSolution") and
        //- fallback dictionary content.
        faSolution
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary& dict
        )
        :
            solution(obr, rOpt, "faSolution", &dict)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("faSolution").
        //  Uses the readOption from the registry.
        explicit faSolution
        (
            const objectRegistry& obr,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, "faSolution", fallback)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("faSolution") and optional contents.
        //  Uses the readOption from the registry.
        faSolution(const objectRegistry& obr, const dictionary& dict)
        :
            solution(obr, "faSolution", &dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

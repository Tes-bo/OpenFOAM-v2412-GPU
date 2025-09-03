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
    Foam::fvSchemes

Description
    Selector class for finite volume differencing schemes.
    fvMesh is derived from fvSchemes so that all fields have access to the
    fvSchemes from the mesh reference they hold.

    The optional fallback dictionary content for constructors is used
    when a file is missing or for a NO_READ, with a null pointer being
    treated like an empty dictionary.

SeeAlso
    Foam::schemesLookup

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvSchemes_H
#define Foam_fvSchemes_H

#include "schemesLookup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvSchemes Declaration
\*---------------------------------------------------------------------------*/

class fvSchemes
:
    public schemesLookup
{
public:

    //- No copy construct
    fvSchemes(const fvSchemes&) = delete;

    //- No copy assignment
    void operator=(const fvSchemes&) = delete;


    // Constructors

        //- Construct for objectRegistry, readOption, (system) dictionary name.
        fvSchemes
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, rOpt, dictName, fallback)
        {}

        //- Construct for objectRegistry, (system) dictionary name
        //- using the readOption from the registry.
        fvSchemes
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, dictName, fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("fvSchemes").
        fvSchemes
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, rOpt, "fvSchemes", fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("fvSchemes") and
        //- fallback dictionary content.
        fvSchemes
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary& dict
        )
        :
            schemesLookup(obr, rOpt, "fvSchemes", &dict)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("fvSchemes").
        //  Uses the readOption from the registry.
        explicit fvSchemes
        (
            const objectRegistry& obr,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, "fvSchemes", fallback)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("fvSchemes") and
        //- fallback dictionary content.
        //  Uses the readOption from the registry.
        fvSchemes(const objectRegistry& obr, const dictionary& dict)
        :
            schemesLookup(obr, "fvSchemes", &dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for fvSchemes - same content for all ranks
template<>
struct is_globalIOobject<fvSchemes> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

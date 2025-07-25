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
    Foam::faSchemes

Description
    Selector class for finite area differencing schemes.
    faMesh is derived from faSchemes so that all fields have access to the
    faSchemes from the mesh reference they hold.

    The optional fallback dictionary content for constructors is used
    when a file is missing or for a NO_READ, with a null pointer being
    treated like an empty dictionary.

SeeAlso
    Foam::schemesLookup

\*---------------------------------------------------------------------------*/

#ifndef Foam_faSchemes_H
#define Foam_faSchemes_H

#include "schemesLookup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faSchemes Declaration
\*---------------------------------------------------------------------------*/

class faSchemes
:
    public schemesLookup
{
public:

    // Generated Methods

        //- No copy construct
        faSchemes(const faSchemes&) = delete;

        //- No copy assignment
        void operator=(const faSchemes&) = delete;


    // Constructors

        //- Construct for objectRegistry, readOption, (system) dictionary name.
        faSchemes
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
        faSchemes
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, dictName, fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("faSchemes").
        faSchemes
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, rOpt, "faSchemes", fallback)
        {}

        //- Construct for objectRegistry, readOption with the
        //- default dictionary name ("faSchemes") and
        //- fallback dictionary content.
        faSchemes
        (
            const objectRegistry& obr,
            IOobjectOption::readOption rOpt,
            const dictionary& dict
        )
        :
            schemesLookup(obr, rOpt, "faSchemes", &dict)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("faSchemes").
        //  Uses the readOption from the registry.
        explicit faSchemes
        (
            const objectRegistry& obr,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, "faSchemes", fallback)
        {}

        //- Construct for objectRegistry with the
        //- default dictionary name ("faSchemes") and
        //- fallback dictionary content.
        //  Uses the readOption from the registry.
        faSchemes(const objectRegistry& obr, const dictionary& dict)
        :
            schemesLookup(obr, "faSchemes", &dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for faSchemes - same content for all ranks
template<>
struct is_globalIOobject<faSchemes> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

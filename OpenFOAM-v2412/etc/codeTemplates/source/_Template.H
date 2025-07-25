/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR,AFFILIATION
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
    Foam::CLASSNAME

Description

SourceFiles
    CLASSNAMEI.H
    CLASSNAME.C
    CLASSNAMEIO.C

\*---------------------------------------------------------------------------*/

#ifndef CLASSNAME_H
#define CLASSNAME_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;
class CLASSNAME;

Istream& operator>>(Istream&, CLASSNAME&);
Ostream& operator<<(Ostream&, const CLASSNAME&);


/*---------------------------------------------------------------------------*\
                         Class CLASSNAME Declaration
\*---------------------------------------------------------------------------*/

class CLASSNAME
:
    public baseClassName
{
    // Private Data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- No copy construct
        CLASSNAME(const CLASSNAME&) = delete;

        //- No copy assignment
        void operator=(const CLASSNAME&) = delete;


public:

    // Static Data Members

        //- Static data staticData
        static const dataType staticData;


    // Generated Methods

//        //- Default construct
//        CLASSNAME() = default;
//
//        //- Copy construct
//        CLASSNAME(const CLASSNAME&) = default;
//
//        //- Copy assignment
//        CLASSNAME& operator=(const CLASSNAME&) = default;


    // Constructors

        //- Default construct
        CLASSNAME();

        //- Construct from components
        CLASSNAME(const dataType& data);

        //- Copy construct
        CLASSNAME(const CLASSNAME&);

        //- Construct from Istream
        explicit CLASSNAME(Istream& is);


    // Selectors

        //- Select default constructed
        static autoPtr<CLASSNAME> New();


    //- Destructor
    ~CLASSNAME();


    // Member Functions

    // Access

    // Check

    // Edit

    // Write


    // Member Operators

        void operator=(const CLASSNAME&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, CLASSNAME&);
        friend Ostream& operator<<(Ostream&, const CLASSNAME&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CLASSNAMEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

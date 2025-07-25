/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::regionProperties

Description
    Simple class to hold region information for coupled region simulations.

    Gives per physics ('fluid', 'solid') the names of the regions. There
    is no assumption on this level that a region should only have a single
    set of physics.

SourceFiles
    regionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionProperties_H
#define Foam_regionProperties_H

#include "HashTable.H"
#include "wordList.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public HashTable<wordList>
{
public:

    // Constructors

        //- Construct from Time
        explicit regionProperties(const Time& runTime);

        //- Construct from Time with specified read option
        regionProperties
        (
            const Time& runTime,
            IOobjectOption::readOption rOpt
        );


    //- Destructor
    ~regionProperties() = default;


    // Member Functions

        //- Total count of all region names.
        label count() const;

        //- The region names. Sorted by region type.
        wordList names() const;

        //- The region names in sorted order.
        wordList sortedNames() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

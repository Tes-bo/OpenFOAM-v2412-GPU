/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 Shell Research Ltd.
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::PDRpatchDef

Description
    Bookkeeping for patch definitions

SourceFiles
    PDRpatchDef.H

\*---------------------------------------------------------------------------*/

#ifndef PDRpatchDef_H
#define PDRpatchDef_H

#include "word.H"
#include "scalar.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class PDRpatchDef Declaration
\*---------------------------------------------------------------------------*/

class PDRpatchDef
{
public:

    //- Patch predefines
    enum predefined
    {
        BLOCKED_FACE = 0,
        MERGING_PATCH = 1,
        WALL_PATCH = 2,
        LAST_PREDEFINED = 2, // First user patch number will be 1 more
        NUM_PREDEFINED = 3
    };

    //- Names for predefined types
    static const Enum<predefined> names;


    // Data Members

    word patchName;

    label patchType;

    scalar blowoffPress;

    scalar blowoffTime;


    // Constructors

        //- Construct null
        PDRpatchDef()
        :
             patchName(),
             patchType(0),
             blowoffPress(0),
             blowoffTime(0)
        {}

        //- Construct with given patch name
        explicit PDRpatchDef(const word& name)
        :
             patchName(name),
             patchType(0),
             blowoffPress(0),
             blowoffTime(0)
        {}


        //- Construct with given patch name
        PDRpatchDef& operator=(const PDRpatchDef&) = default;
        PDRpatchDef& operator=(PDRpatchDef&&) = default;

        //- Assign new patch name
        void operator=(const std::string& newName);
};


typedef PDRpatchDef PATCH;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

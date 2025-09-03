/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::surfZoneIOList

Description
    IOobject for a surfZoneList

SourceFiles
    surfZoneIOList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfZoneIOList_H
#define Foam_surfZoneIOList_H

#include "surfZoneList.H"
#include "regIOobject.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class surfZoneIOList Declaration
\*---------------------------------------------------------------------------*/

class surfZoneIOList
:
    public regIOobject,
    public surfZoneList
{
    // Private Member Functions

        //- Return true if contents were read
        //- (controlled by IOobject readOption flags).
        bool readIOcontents();


public:

    //- Runtime type information
    TypeNameNoDebug("surfZoneList");


    // Constructors

        //- Construct from IOobject
        explicit surfZoneIOList(const IOobject& io);

        //- Construct from IOobject and surfZoneList content
        surfZoneIOList(const IOobject& io, const UList<surfZone>& content);

        //- Construct from IOobject and transferring surfZoneList content
        surfZoneIOList(const IOobject& io, surfZoneList&& content);


    //- Destructor
    virtual ~surfZoneIOList() = default;


    // Member Functions

        //- The writeData function - required by regIOobject
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Copy assignment of entries
        void operator=(const surfZoneIOList& rhs);

        //- Copy or move assignment of entries
        using surfZoneList::operator=;
};


// Global Operators

//- Write surfZone list using its writeData() method
Ostream& operator<<(Ostream& os, const surfZoneIOList& zones);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

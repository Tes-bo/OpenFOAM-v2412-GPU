/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::searchableSurfaceModifier::autoPatch

Description
    Patchify triangles based on a feature angle. Replacement of surfaceAutoPatch

SourceFiles
    autoPatch.C

\*---------------------------------------------------------------------------*/

#ifndef autoPatch_H
#define autoPatch_H

#include "searchableSurfaceModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace searchableSurfaceModifiers
{

/*---------------------------------------------------------------------------*\
                          Class autoPatch Declaration
\*---------------------------------------------------------------------------*/

class autoPatch
:
    public searchableSurfaceModifier
{
    // Private data

        //- Feature angle
        const scalar featureAngle_;


public:

    //- Runtime type information
    TypeName("autoPatch");


    // Constructors

        //- Construct from dictionary
        autoPatch(const searchableSurfaces& geometry, const dictionary& dict);

        //- Clone
        autoPtr<searchableSurfaceModifier> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    //- Destructor
    virtual ~autoPatch() = default;


    // Member Functions

        //- Apply this selector
        virtual bool modify(const labelList& regions, searchableSurface&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace searchableSurfaceModifiers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::mixerFvMesh

Description
    A rotating slider mesh

SourceFiles
    mixerFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mixerFvMesh_H
#define mixerFvMesh_H

#include "topoChangerFvMesh.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixerFvMesh Declaration
\*---------------------------------------------------------------------------*/

class mixerFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

        //- Cylindrical coordinate system
        coordSystem::cylindrical csys_;

        // - Rotational speed in revolutions per minute (rpm)
        scalar rpm_;

        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;


    // Private Member Functions

        //- No copy construct
        mixerFvMesh(const mixerFvMesh&) = delete;

        //- No copy assignment
        void operator=(const mixerFvMesh&) = delete;


        //- Add mixer zones and modifiers
        void addZonesAndModifiers();

        //- Calculate moving masks
        void calcMovingMasks() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;


public:

    //- Runtime type information
    TypeName("mixerFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit mixerFvMesh(const IOobject& io);


    //- Destructor
    virtual ~mixerFvMesh();


    // Member Functions

        //- Return coordinate system
        const coordSystem::cylindrical& cs() const
        {
            return csys_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

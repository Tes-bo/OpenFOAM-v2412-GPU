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
    Foam::emptyPointPatch

Description
    Empty plane patch.

SourceFiles
    emptyPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef emptyPointPatch_H
#define emptyPointPatch_H

#include "facePointPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class emptyPointPatch Declaration
\*---------------------------------------------------------------------------*/

class emptyPointPatch
:
    public facePointPatch
{

public:

    //- Runtime type information
    TypeName(emptyPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        emptyPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            facePointPatch(patch, bm)
        {}

        //- Construct given the original patch and a map
        emptyPointPatch
        (
            const emptyPointPatch& patch,
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        )
        :
            facePointPatch(patch, bm, index, mapAddressing, reversePointMap)
        {}

        //- Construct and return a subset clone,
        //- resetting the point list and boundary mesh
        virtual autoPtr<pointPatch> clone
        (
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        ) const
        {
            return autoPtr<pointPatch>::NewFrom<emptyPointPatch>
            (
                *this,
                bm,
                index,
                mapAddressing,
                reversePointMap
            );
        }


    // Member Functions

        //- Return the constraint type this pointPatch implements.
        virtual const word& constraintType() const
        {
            return type();
        }

        //- Accumulate the effect of constraint direction of this patch
        virtual void applyConstraint
        (
            const label pointi,
            pointConstraint&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

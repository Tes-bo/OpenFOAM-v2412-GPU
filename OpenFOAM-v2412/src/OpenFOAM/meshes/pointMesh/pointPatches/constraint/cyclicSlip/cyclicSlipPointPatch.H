/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::cyclicSlipPointPatch

Description
    Cyclic patch with slip constraint

SourceFiles
    cyclicSlipPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicSlipPointPatch_H
#define cyclicSlipPointPatch_H

#include "cyclicPointPatch.H"
#include "cyclicSlipPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class cyclicSlipPointPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicSlipPointPatch
:
    public cyclicPointPatch
{

public:

    //- Runtime type information
    TypeName(cyclicSlipPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        cyclicSlipPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            cyclicPointPatch(patch, bm)
        {}

        //- Construct given the original patch and a map
        cyclicSlipPointPatch
        (
            const cyclicSlipPointPatch& patch,
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        )
        :
            cyclicPointPatch(patch, bm, index, mapAddressing, reversePointMap)
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
            return autoPtr<pointPatch>::NewFrom<cyclicSlipPointPatch>
            (
                *this,
                bm,
                index,
                mapAddressing,
                reversePointMap
            );
        }


    //- Destructor
    virtual ~cyclicSlipPointPatch() = default;


    // Member Functions

        //- Return point unit normals.
        virtual const vectorField& pointNormals() const;

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

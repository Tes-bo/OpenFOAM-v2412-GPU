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
    Foam::coupledFacePointPatch

Description
    coupled patch for post-processing.  Used as the base class for processor
    and cyclic pointPatches

SourceFiles
    coupledFacePointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_coupledFacePointPatch_H
#define Foam_coupledFacePointPatch_H

#include "coupledPointPatch.H"
#include "facePointPatch.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledFacePointPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFacePointPatch
:
    public facePointPatch,
    public coupledPointPatch
{
    // Private Data

        const coupledPolyPatch& coupledPolyPatch_;


protected:

    // Protected Member Functions

        //- Calculate mesh points
        virtual void calcGeometry(PstreamBuffers&) = 0;


        //- No copy construct
        coupledFacePointPatch(const coupledFacePointPatch&) = delete;

        //- No copy assignment
        void operator=(const coupledFacePointPatch&) = delete;


public:

    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledFacePointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );

        //- Construct given the original patch and a map
        coupledFacePointPatch
        (
            const coupledFacePointPatch& pp,
            const pointBoundaryMesh& bm,
            const label index,
            const labelUList& mapAddressing,
            const labelUList& reversePointMap
        );


    //- Destructor
    virtual ~coupledFacePointPatch() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

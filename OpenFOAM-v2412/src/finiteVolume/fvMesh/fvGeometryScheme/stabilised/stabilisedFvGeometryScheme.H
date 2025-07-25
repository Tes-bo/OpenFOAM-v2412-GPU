/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::stabilisedFvGeometryScheme

Description
    Geometry calculation scheme that implements face geometry calculation
    using normal-component-of-area weighted triangle contributions.

    This implements the Foundation 'Corrected face-centre calculations'
    as a separate geometry scheme. Only implements the primitiveMesh parts,
    not the individual face calculation.

SourceFiles
    stabilisedFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef stabilisedFvGeometryScheme_H
#define stabilisedFvGeometryScheme_H

#include "basicFvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
               Class stabilisedFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class stabilisedFvGeometryScheme
:
    public basicFvGeometryScheme
{
protected:

    // Protected Member Functions

        //- Calculate face area and centre weighted using pyramid height
        static void makeFaceCentresAndAreas
        (
            const polyMesh& mesh,
            const pointField& p,
            vectorField& fCtrs,
            vectorField& fAreas
        );


private:

    // Private Member Functions

        //- No copy construct
        stabilisedFvGeometryScheme
        (
            const stabilisedFvGeometryScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const stabilisedFvGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("stabilised");


    // Constructors

        //- Construct from mesh
        stabilisedFvGeometryScheme
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~stabilisedFvGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();

        //- Calculate geometry quantities using mesh topology and provided
        //- points. If oldPoints provided only does local update. Returns
        //- true if anything changed, false otherwise
        virtual bool updateGeom
        (
            const pointField& points,
            const refPtr<pointField>& oldPoints,    // optional old points
            pointField& faceCentres,
            vectorField& faceAreas,
            pointField& cellCentres,
            scalarField& cellVolumes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

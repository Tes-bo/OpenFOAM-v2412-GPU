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
    Foam::basicFvGeometryScheme

Description
    Default geometry calculation scheme. Slight stabilisation for bad meshes.

SourceFiles
    basicFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef basicFvGeometryScheme_H
#define basicFvGeometryScheme_H

#include "fvGeometryScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class basicFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class basicFvGeometryScheme
:
    public fvGeometryScheme
{
    // Private Member Functions

        //- No copy construct
        basicFvGeometryScheme(const basicFvGeometryScheme&) = delete;

        //- No copy assignment
        void operator=(const basicFvGeometryScheme&) = delete;


public:

    //- Runtime type information
    TypeName("basic");


    // Constructors

        //- Construct from mesh
        basicFvGeometryScheme(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~basicFvGeometryScheme() = default;


    // Member Functions

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();

        //- Return linear difference weighting factors
        virtual tmp<surfaceScalarField> weights() const;

        //- Return cell-centre difference coefficients
        virtual tmp<surfaceScalarField> deltaCoeffs() const;

        //- Return non-orthogonal cell-centre difference coefficients
        virtual tmp<surfaceScalarField> nonOrthDeltaCoeffs() const;

        //- Return non-orthogonality correction vectors
        virtual tmp<surfaceVectorField> nonOrthCorrectionVectors() const;

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

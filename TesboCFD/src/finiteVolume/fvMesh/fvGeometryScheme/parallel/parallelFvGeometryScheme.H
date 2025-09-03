/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::parallelFvGeometryScheme

Description
    Geometry calculation scheme with explicit sync of face-geometry
    across processor patches.

    - applies face centre and normal (negated) to the neighbour side
    of a processor boundary.
    - recalculates the resulting cell volumes.
    - mainly interesting in single-precision in that it removes the
    different truncation error from circulating in different order. This
    can cause problems when calculating global transformations.
    - it does change the face normal/centre so might affect the
    cell-closedness of the attached cells.

Usage
    In system/fvSchemes add a geometry scheme:
    \verbatim
    geometry
    {
        type            parallel;

        // Optional underlying geometry scheme. (default is basic)
        geometry
        {
            type        solidBody;
        }
    }
    \endverbatim

    The default behaviour will print the total number of adjusted face and
    cell properties:
    \verbatim
        parallelFvGeometryScheme::movePoints() :
            adjusted geometry of faces:113 of cells:54
    \endverbatim
    In addition a warning will be printed if an internal face gives a different
    result.

SourceFiles
    parallelFvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_parallelFvGeometryScheme_H
#define Foam_parallelFvGeometryScheme_H

#include "basicFvGeometryScheme.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class parallelFvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class parallelFvGeometryScheme
:
    public fvGeometryScheme
{
protected:

    // Protected Data

        //- Dictionary for underlying scheme
        dictionary dict_;

        //- Demand-driven construction of underlying scheme
        mutable tmp<fvGeometryScheme> geometryPtr_;


private:

    // Private Member Functions

        //- Swap processor-face geometry
        void adjustGeometry
        (
            pointField& faceCentres,
            vectorField& faceAreas
        ) const;

        //- Swap processor-face geometry
        void adjustGeometry();

        //- No copy construct
        parallelFvGeometryScheme(const parallelFvGeometryScheme&) = delete;

        //- No copy assignment
        void operator=(const parallelFvGeometryScheme&) = delete;


public:

    // Public classes

        //- Transformation for oriented fields
        class transformOriented
        {
        public:

            template<class Type>
            void operator()(const coupledPolyPatch& cpp, UList<Type>& fld) const
            {
                if (cpp.parallel())
                {
                    for (auto& f : fld)
                    {
                        f = -f;
                    }
                }
                else
                {
                    transformList(-cpp.forwardT(), fld);
                }
            }
        };


public:

    //- Runtime type information
    TypeName("parallel");


    // Constructors

        //- Construct from mesh
        parallelFvGeometryScheme
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~parallelFvGeometryScheme() = default;


    // Member Functions

        //- Construct underlying fvGeometryScheme
        const fvGeometryScheme& geometry() const;

        //- Do what is necessary if the mesh has moved
        virtual void movePoints();

        //- Update mesh for topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);

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

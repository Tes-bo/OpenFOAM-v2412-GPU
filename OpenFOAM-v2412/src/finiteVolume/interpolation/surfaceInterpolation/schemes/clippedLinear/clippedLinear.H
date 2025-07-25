/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::clippedLinear

Group
    grpFvSurfaceInterpolationSchemes

Description
    Central-differencing interpolation scheme using clipped-weights to
    improve stability on meshes with very rapid variations in cell size.

SourceFiles
    clippedLinear.C

\*---------------------------------------------------------------------------*/

#ifndef clippedLinear_H
#define clippedLinear_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class clippedLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class clippedLinear
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        scalar wfLimit_;


    // Private Member Functions

        static scalar calcWfLimit(scalar cellSizeRatio)
        {
            if (cellSizeRatio <= 0 || cellSizeRatio > 1)
            {
                FatalErrorInFunction
                    << "Given cellSizeRatio of " << cellSizeRatio
                    << " is not between 0 and 1"
                    << exit(FatalError);
            }

            return cellSizeRatio/(1.0 + cellSizeRatio);
        }


        //- No copy assignment
        void operator=(const clippedLinear&) = delete;


public:

    //- Runtime type information
    TypeName("clippedLinear");


    // Constructors

        //- Construct from mesh and cellSizeRatio
        clippedLinear(const fvMesh& mesh, const scalar cellSizeRatio)
        :
            surfaceInterpolationScheme<Type>(mesh),
            wfLimit_(calcWfLimit(cellSizeRatio))
        {}

        //- Construct from Istream
        clippedLinear(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            wfLimit_(calcWfLimit(readScalar(is)))
        {}

        //- Construct from faceFlux and Istream
        clippedLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            wfLimit_(calcWfLimit(readScalar(is)))
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            const fvMesh& mesh = this->mesh();

            tmp<surfaceScalarField> tcdWeights
            (
                mesh.surfaceInterpolation::weights()
            );
            const surfaceScalarField& cdWeights = tcdWeights();

            const scalarMinMax bounds(wfLimit_, 1 - wfLimit_);

            tmp<surfaceScalarField> tclippedLinearWeights
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "clippedLinearWeights",
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimless
                )
            );
            surfaceScalarField& clippedLinearWeights =
                tclippedLinearWeights.ref();

            clippedLinearWeights.primitiveFieldRef() =
                clamp(cdWeights.primitiveField(), bounds);

            surfaceScalarField::Boundary& clwbf =
                clippedLinearWeights.boundaryFieldRef();

            forAll(mesh.boundary(), patchi)
            {
                if (clwbf[patchi].coupled())
                {
                    clwbf[patchi] =
                        clamp(cdWeights.boundaryField()[patchi], bounds);
                }
                else
                {
                    clwbf[patchi] = cdWeights.boundaryField()[patchi];
                }
            }

            return tclippedLinearWeights;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::blendedEdgeInterpolation

Description
    linear/upwind blended differencing scheme.

SourceFiles
    blendedEdgeInterpolationMake.C

\*---------------------------------------------------------------------------*/

#ifndef blendedEdgeInterpolation_H
#define blendedEdgeInterpolation_H

#include "linearEdgeInterpolation.H"
#include "upwindEdgeInterpolation.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class blendedEdgeInterpolation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class blendedEdgeInterpolation
:
    public linearEdgeInterpolation<Type>,
    public upwindEdgeInterpolation<Type>
{
    // Private Data

        //- Blending factor
        const scalar blendingFactor_;


public:

    //- Runtime type information
    TypeName("blended");


    // Generated Methods

        //- No copy construct
        blendedEdgeInterpolation(const blendedEdgeInterpolation&) = delete;

        //- No copy assignment
        void operator=(const blendedEdgeInterpolation&) = delete;


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        blendedEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            const scalar blendingFactor
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh),
            upwindEdgeInterpolation<Type>(mesh, faceFlux),
            blendingFactor_(blendingFactor)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the database
        blendedEdgeInterpolation
        (
            const faMesh& mesh,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh),
            upwindEdgeInterpolation<Type>
            (
                mesh,
                mesh.thisDb().lookupObject<edgeScalarField>
                (
                    word(is)
                )
            ),
            blendingFactor_(readScalar(is))
        {}

        //- Construct from mesh, faceFlux and Istream
        blendedEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            linearEdgeInterpolation<Type>(mesh),
            upwindEdgeInterpolation<Type>(mesh, faceFlux),
            blendingFactor_(readScalar(is))
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>& vf
        ) const
        {
            return
                blendingFactor_*
                linearEdgeInterpolation<Type>::weights(vf)
              + (1 - blendingFactor_)*
                upwindEdgeInterpolation<Type>::weights(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

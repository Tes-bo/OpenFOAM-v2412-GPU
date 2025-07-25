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
    Foam::upwindEdgeInterpolation

Description
    Upwind differencing scheme class.

SourceFiles
    upwindEdgeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef upwindEdgeInterpolation_H
#define upwindEdgeInterpolation_H

#include "edgeInterpolationScheme.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class upwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class upwindEdgeInterpolation
:
    virtual public edgeInterpolationScheme<Type>
{
    // Private Data

        //- Face flux
        const edgeScalarField& faceFlux_;


public:

    //- Runtime type information
    TypeName("upwind");


    // Generated Methods

        //- No copy construct
        upwindEdgeInterpolation(const upwindEdgeInterpolation&) = delete;

        //- No copy assignment
        void operator=(const upwindEdgeInterpolation&) = delete;


    // Constructors

        //- Construct from faceFlux
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}

        //- Construct from Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the database
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            Istream& is
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_
            (
                mesh.thisDb().lookupObject<edgeScalarField>
                (
                    word(is)
                )
            )
        {}

        //- Construct from faceFlux and Istream
        upwindEdgeInterpolation
        (
            const faMesh& mesh,
            const edgeScalarField& faceFlux,
            Istream&
        )
        :
            edgeInterpolationScheme<Type>(mesh),
            faceFlux_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<edgeScalarField> weights
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const
        {
            return pos(faceFlux_);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::fa::gaussGrad

Description
    Basic second-order gradient scheme using face-interpolation
    and Gauss' theorem.

SourceFiles
    gaussFaGrad.C
    gaussFaGrads.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gaussFaGrad_H
#define Foam_gaussFaGrad_H

#include "faGradScheme.H"
#include "edgeInterpolationScheme.H"
#include "linearEdgeInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                          Class gaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussGrad
:
    public fa::gradScheme<Type>
{
    // Private Data

        //- Interpolation scheme
        tmp<edgeInterpolationScheme<Type>> tinterpScheme_;


public:

    //- Runtime type information
    TypeName("Gauss");


    // Generated Methods

        //- No copy construct
        gaussGrad(const gaussGrad&) = delete;

        //- No copy assignment
        void operator=(const gaussGrad&) = delete;


    // Constructors

        //- Construct from mesh
        gaussGrad(const faMesh& mesh)
        :
            gradScheme<Type>(mesh),
            tinterpScheme_(new linearEdgeInterpolation<Type>(mesh))
        {}

        //- Construct from Istream
        gaussGrad(const faMesh& mesh, Istream& is)
        :
            gradScheme<Type>(mesh)
        {
            if (is.eof())
            {
                tinterpScheme_.reset
                (
                    new linearEdgeInterpolation<Type>(mesh)
                );
            }
            else
            {
                tinterpScheme_.reset
                (
                    edgeInterpolationScheme<Type>::New(mesh, is)
                );
            }
        }


    // Member Functions

        //- Return the gradient of the given field
        //- calculated using Gauss' theorem on the given surface field
        static
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > gradf
        (
            const GeometricField<Type, faePatchField, edgeMesh>&,
            const word& name
        );

        //- Return the gradient of the given field to the
        //- gradScheme::grad for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > calcGrad
        (
            const GeometricField<Type, faPatchField, areaMesh>& vsf,
            const word& name
        ) const;


        //- Correct the boundary values of the gradient using the patchField
        // snGrad functions
        static void correctBoundaryConditions
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussFaGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::fv::iterativeGaussGrad

Description
    A second-order gradient scheme using face-interpolation,
    Gauss' theorem and iterative skew correction.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    gradSchemes
    {
        grad(<term>)   iterativeGauss <interpolation scheme> <number of iters>;
    }
    \endverbatim

    and by using \c system/fvSolution:
    \verbatim
    relaxationFactors
    {
        fields
        {
            grad(<term>)    <relaxation factor>;
        }
    }
    \endverbatim

SourceFiles
    iterativeGaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef iterativeGaussGrad_H
#define iterativeGaussGrad_H

#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class iterativeGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class iterativeGaussGrad
:
    public fv::gaussGrad<Type>
{
    // Private Data

        //- Number of skew-correction iterations
        label nIter_;


    // Private Member Functions

        //- No copy construct
        iterativeGaussGrad(const iterativeGaussGrad&) = delete;

        //- No copy assignment
        void operator=(const iterativeGaussGrad&) = delete;


public:

    //- Runtime type information
    TypeName("iterativeGauss");


    // Constructors

        //- Construct from mesh
        iterativeGaussGrad(const fvMesh& mesh)
        :
            gaussGrad<Type>(mesh),
            nIter_(1)
        {}

        //- Construct from mesh and Istream
        iterativeGaussGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gaussGrad<Type>(mesh, schemeData),
            nIter_(readLabel(schemeData))
        {
            if (nIter_ <= 0)
            {
                FatalIOErrorInFunction(schemeData)
                    << "nIter = " << nIter_
                    << " should be > 0"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        //- Return the gradient of the given field
        //- to the gradScheme::grad for optional caching
        virtual tmp
        <
            GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "iterativeGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

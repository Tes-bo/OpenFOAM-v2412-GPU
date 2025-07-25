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
    Foam::fa::gaussLaplacianScheme

Description
    Basic second-order laplacian using face-gradients and Gauss' theorem.

SourceFiles
    gaussFaLaplacianScheme.C
    gaussFaLaplacianSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_gaussFaLaplacianScheme_H
#define Foam_gaussFaLaplacianScheme_H

#include "faLaplacianScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                    Class gaussLaplacianScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gaussLaplacianScheme
:
    public fa::laplacianScheme<Type>
{
public:

    //- Runtime type information
    TypeName("Gauss");


    // Generated Methods

        //- No copy construct
        gaussLaplacianScheme(const gaussLaplacianScheme&) = delete;

        //- No copy assignment
        void operator=(const gaussLaplacianScheme&) = delete;


    // Constructors

        //- Construct from mesh
        gaussLaplacianScheme(const faMesh& mesh)
        :
            laplacianScheme<Type>(mesh)
        {}

        //- Construct from Istream
        gaussLaplacianScheme(const faMesh& mesh, Istream& is)
        :
            laplacianScheme<Type>(mesh, is)
        {}


    // Member Functions

        tmp<faMatrix<Type>> famLaplacian
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facLaplacian
        (
            const edgeScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gaussFaLaplacianScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::fa::steadyStateDdtScheme

Description
    SteadyState implicit/explicit ddt which returns 0.

SourceFiles
    steadyStateFaDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef steadyStateFaDdtScheme_H
#define steadyStateFaDdtScheme_H

#include "faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                   Class steadyStateFaDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class steadyStateFaDdtScheme
:
    public fa::faDdtScheme<Type>
{
public:

    //- Runtime type information
    TypeName("steadyState");


    // Generated Methods

        //- No copy construct
        steadyStateFaDdtScheme(const steadyStateFaDdtScheme&) = delete;

        //- No copy assignment
        void operator=(const steadyStateFaDdtScheme&) = delete;


    // Constructors

        //- Construct from mesh
        steadyStateFaDdtScheme(const faMesh& mesh)
        :
            faDdtScheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        steadyStateFaDdtScheme(const faMesh& mesh, Istream& is)
        :
            faDdtScheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return fa::faDdtScheme<Type>::mesh();
        }

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "steadyStateFaDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

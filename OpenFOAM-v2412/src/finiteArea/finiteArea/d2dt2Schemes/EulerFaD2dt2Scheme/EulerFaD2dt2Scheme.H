/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 Volkswagen AG
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
    Foam::fa::EulerFaD2dt2Scheme

Description
    First-order Euler implicit d2dt2 using the current and two previous
    time-step values.

SourceFiles
    EulerFaD2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef EulerFaD2dt2Scheme_H
#define EulerFaD2dt2Scheme_H

#include "faD2dt2Scheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                      Class EulerD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class EulerFaD2dt2Scheme
:
    public fa::faD2dt2Scheme<Type>
{
    // Private Member Functions

        //- Return the current time-step
        scalar deltaT_() const;

        //- Return the previous time-step
        scalar deltaT0_() const;


public:

    //- Runtime type information
    TypeName("Euler");


    // Generated Methods

        //- No copy construct
        EulerFaD2dt2Scheme(const EulerFaD2dt2Scheme&) = delete;

        //- No copy assignment
        void operator=(const EulerFaD2dt2Scheme&) = delete;


    // Constructors

        //- Construct from mesh
        EulerFaD2dt2Scheme(const faMesh& mesh)
        :
            faD2dt2Scheme<Type>(mesh)
        {}

        //- Construct from mesh and Istream
        EulerFaD2dt2Scheme(const faMesh& mesh, Istream& is)
        :
            faD2dt2Scheme<Type>(mesh, is)
        {}


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const
        {
            return fa::faD2dt2Scheme<Type>::mesh();
        }

        tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const dimensioned<Type>
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famD2dt2
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        tmp<faMatrix<Type>> famD2dt2
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
#   include "EulerFaD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

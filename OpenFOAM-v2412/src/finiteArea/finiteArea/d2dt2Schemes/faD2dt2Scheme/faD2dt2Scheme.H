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
    Foam::fa::faD2dt2Scheme

Description
    Abstract base class for d2dt2 schemes.

SourceFiles
    faD2dt2Scheme.C

\*---------------------------------------------------------------------------*/

#ifndef faD2dt2Scheme_H
#define faD2dt2Scheme_H

#include "dimensionedType.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Type>
class faMatrix;

class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class faD2dt2Scheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faD2dt2Scheme
:
    public refCount
{
protected:

    // Protected Data

        //- Reference to mesh
        const faMesh& mesh_;


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faD2dt2Scheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Generated Methods

        //- No copy construct
        faD2dt2Scheme(const faD2dt2Scheme&) = delete;

        //- No copy assignment
        void operator=(const faD2dt2Scheme&) = delete;


    // Constructors

        //- Construct from mesh
        faD2dt2Scheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from mesh and Istream
        faD2dt2Scheme(const faMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new d2dt2Scheme created on freestore
        static tmp<faD2dt2Scheme<Type>> New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~faD2dt2Scheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const noexcept { return mesh_; }

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const dimensioned<Type>
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facD2dt2
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famD2dt2
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famD2dt2
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famD2dt2
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaD2dt2TypeScheme(SS, Type)                                        \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::fa::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fa                                                           \
        {                                                                      \
            faD2dt2Scheme<Type>::addIstreamConstructorToTable<SS<Type>>        \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }


#define makeFaD2dt2Scheme(SS)                                                  \
                                                                               \
makeFaD2dt2TypeScheme(SS, scalar)                                              \
makeFaD2dt2TypeScheme(SS, vector)                                              \
makeFaD2dt2TypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faD2dt2Scheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

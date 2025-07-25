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
    Foam::fa::DdtScheme

Description
    Abstract base class for finite volume calculus ddt schemes.

SourceFiles
    faDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef faDdtScheme_H
#define faDdtScheme_H

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
                           Class faDdtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faDdtScheme
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
            faDdtScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Generated Methods

        //- No copy construct
        faDdtScheme(const faDdtScheme&) = delete;

        //- No copy assignment
        void operator=(const faDdtScheme&) = delete;


    // Constructors

        //- Construct from mesh
        faDdtScheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from mesh and Istream
        faDdtScheme(const faMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new faDdtScheme created on freestore
        static tmp<faDdtScheme<Type>> New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~faDdtScheme();


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const noexcept { return mesh_; }

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const dimensioned<Type>
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const dimensioned<Type>
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<GeometricField<Type, faPatchField, areaMesh>> facDdt0
        (
            const areaScalarField&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famDdt
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famDdt
        (
            const dimensionedScalar&,
            const GeometricField<Type, faPatchField, areaMesh>&
        ) = 0;

        virtual tmp<faMatrix<Type>> famDdt
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

#define makeFaDdtTypeScheme(SS, Type)                                          \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::fa::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fa                                                           \
        {                                                                      \
            faDdtScheme<Type>::addIstreamConstructorToTable<SS<Type>>          \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }


#define makeFaDdtScheme(SS)                                                    \
                                                                               \
makeFaDdtTypeScheme(SS, scalar)                                                \
makeFaDdtTypeScheme(SS, vector)                                                \
makeFaDdtTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

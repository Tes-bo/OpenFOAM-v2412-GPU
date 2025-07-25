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
    Foam::fa::gradScheme

Description
    Abstract base class for finite area calculus gradient schemes.

SourceFiles
    faGradScheme.C
    faGradSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef faGradScheme_H
#define faGradScheme_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                         Class gradScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class gradScheme
:
    public refCount
{
    // Private Data

        //- Reference to mesh
        const faMesh& mesh_;


public:

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            gradScheme,
            Istream,
            (const faMesh& mesh, Istream& schemeData),
            (mesh, schemeData)
        );


    // Generated Methods

        //- No copy construct
        gradScheme(const gradScheme&) = delete;

        //- No copy assignment
        void operator=(const gradScheme&) = delete;


    // Constructors

        //- Construct from mesh
        gradScheme(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return a pointer to a new gradScheme created on freestore
        static tmp<gradScheme<Type>> New
        (
            const faMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~gradScheme() = default;


    // Member Functions

        //- Return mesh reference
        const faMesh& mesh() const noexcept { return mesh_; }

        //- Calculate and return the grad of the given field.
        //  Used by grad either to recalculate the cached gradient when it is
        //  out of date with respect to the field or when it is not cached.
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > calcGrad
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const word& name
        ) const = 0;

        //- Calculate and return the grad of the given field
        //- which may have been cached
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            const word& name
        ) const;

        //- Calculate and return the grad of the given field
        //- with the default name
        //- which may have been cached
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        ) const;

        //- Calculate and return the grad of the given field
        //- with the default name
        //- which may have been cached
        tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > grad
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFaGradTypeScheme(SS, Type)                                         \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::fa::SS<Foam::Type>, 0);          \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace fa                                                           \
        {                                                                      \
            gradScheme<Type>::addIstreamConstructorToTable<SS<Type>>          \
                add##SS##Type##IstreamConstructorToTable_;                     \
        }                                                                      \
    }


#define makeFaGradScheme(SS)                                                   \
                                                                               \
makeFaGradTypeScheme(SS, scalar)                                               \
makeFaGradTypeScheme(SS, vector)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::PatchFunction1

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or time varying, and
    provide functions to return the (interpolated) value, and integral between
    limits.

    Extends the Function1 class by adding autoMap and rMap functions

SourceFiles
    PatchFunction1.C
    PatchFunction1New.C

SeeAlso
    Foam::Function1

\*---------------------------------------------------------------------------*/

#ifndef Foam_PatchFunction1_H
#define Foam_PatchFunction1_H

#include "patchFunction1Base.H"
#include "coordinateScaling.H"
#include "Field.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
template<class Type> class PatchFunction1;

template<class Type>
Ostream& operator<<(Ostream&, const PatchFunction1<Type>&);

/*---------------------------------------------------------------------------*\
                       Class PatchFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class PatchFunction1
:
    public patchFunction1Base
{
    // Private Member Functions

        //- Selector, with alternative entry etc
        static autoPtr<PatchFunction1<Type>> New
        (
            const polyPatch& pp,
            const word& entryName,
            const entry* eptr,
            const dictionary& dict,
            const bool faceValues,
            const bool mandatory
        );


protected:

    // Protected Data

        //- Optional local coordinate system and scaling
        coordinateScaling<Type> coordSys_;


    // Protected Member Functions

        //- No copy assignment
        void operator=(const PatchFunction1<Type>&) = delete;


public:

    // Data Types

        //- The return type is a field of values
        typedef Field<Type> returnType;


    //- Runtime type information
    TypeName("PatchFunction1")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        PatchFunction1,
        dictionary,
        (
            const polyPatch& pp,
            const word& type,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues
        ),
        (pp, type, entryName, dict, faceValues)
    );


    // Constructors

        //- Construct from polyPatch and entry name
        PatchFunction1
        (
            const polyPatch& pp,
            const word& entryName,
            const bool faceValues = true
        );

        //- Construct from polyPatch, dictionary and entry name
        PatchFunction1
        (
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct setting patch
        explicit PatchFunction1
        (
            const PatchFunction1<Type>& rhs,
            const polyPatch& pp
        );

        //- Copy construct
        explicit PatchFunction1(const PatchFunction1<Type>& rhs);

        //- Return a clone
        virtual tmp<PatchFunction1<Type>> clone() const = 0;

        //- Return a clone, setting patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const = 0;


    // Factory Methods

        //- Clone a PatchFunction1
        template<class Derived>
        static tmp<PatchFunction1<Type>>
        Clone(const Derived& fun)
        {
            return tmp<PatchFunction1<Type>>(new Derived(fun));
        }

        //- Clone a PatchFunction1 with a patch
        template<class Derived>
        static tmp<PatchFunction1<Type>>
        Clone(const Derived& fun, const polyPatch& pp)
        {
            return tmp<PatchFunction1<Type>>(new Derived(fun, pp));
        }


        //- Selector
        static autoPtr<PatchFunction1<Type>> New
        (
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true,
            const bool mandatory = true
        );

        //- Compatibility selector
        static autoPtr<PatchFunction1<Type>> NewCompat
        (
            const polyPatch& pp,
            const word& entryName,
            std::initializer_list<std::pair<const char*,int>> compat,
            const dictionary& dict,
            const bool faceValues = true,
            const bool mandatory = true
        );

        //- An optional selector
        static autoPtr<PatchFunction1<Type>> NewIfPresent
        (
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );


    // Caching Selectors - accept wildcards in dictionary

        //- Selector with external storage.
        //- This also allows wildcard matches in a dictionary
        static refPtr<PatchFunction1<Type>> New
        (
            HashPtrTable<PatchFunction1<Type>>& cache,
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            enum keyType::option matchOpt = keyType::LITERAL,
            const bool faceValues = true,
            const bool mandatory = true
        );


    //- Destructor
    virtual ~PatchFunction1() = default;


    // Member Functions

        //- Is value constant (i.e. independent of x)
        virtual bool constant() const { return false; }

        //- Can function be evaluated?
        virtual bool good() const { return true; }

        //- Is value uniform (i.e. independent of coordinate)
        virtual bool uniform() const = 0;


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual tmp<Field<Type>> value(const scalar x) const;

            //- Integrate between two (scalar) values
            virtual tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;

            //- Helper: optionally convert coordinates to local coordinates
            virtual tmp<pointField> localPosition
            (
                const pointField& globalPos
            ) const;

            //- Apply optional transformation
            virtual tmp<Field<Type>> transform(const Field<Type>& fld) const;

            //- Apply optional transformation
            virtual tmp<Field<Type>> transform
            (
                const tmp<Field<Type>>& tfld
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<Type>& rhs,
                const labelList& addr
            );


        // I/O

            //- Ostream Operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const PatchFunction1<Type>& rhs
            );

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define PatchFunction1 run-time selection
#define makePatchFunction1(Type)                                               \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(PatchFunction1<Type>, 0);              \
                                                                               \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        PatchFunction1<Type>,                                                  \
        dictionary                                                             \
    );


// Define (templated) PatchFunction1, add to (templated) run-time selection
#define makePatchFunction1Type(SS, Type)                                       \
                                                                               \
    defineNamedTemplateTypeNameAndDebug(PatchFunction1Types::SS<Type>, 0);     \
                                                                               \
    PatchFunction1<Type>::adddictionaryConstructorToTable                      \
        <PatchFunction1Types::SS<Type>>                                        \
        add##SS##Type##ConstructorToTable_;


// Define (non-templated) PatchFunction1, add to (templated) run-time selection
#define makeConcretePatchFunction1Type(SS, Type)                               \
                                                                               \
    defineTypeNameAndDebug(SS, 0);                                             \
                                                                               \
    PatchFunction1<Type>::adddictionaryConstructorToTable                      \
        <PatchFunction1Types::SS>                                              \
        add##SS##Type##ConstructorToTable_;


// Define scalar PatchFunction1 and add to (templated) run-time selection
#define makeScalarPatchFunction1(SS)                                           \
                                                                               \
    makeConcretePatchFunction1Type(SS, scalar);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PatchFunction1.C"
    #include "PatchFunction1New.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::PatchFunction1Types::ConstantField

Description
    Templated function that returns a constant value.

    Usage - for entry \<entryName\> returning the value \<value\>:
    \verbatim
        <entryName>    constant  <value>;
    \endverbatim
    or
    \verbatim
        <entryName>
        {
            type    constant;
            value   <value>;
        }
    \endverbatim

SourceFiles
    ConstantField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PatchFunction1Types_ConstantField_H
#define Foam_PatchFunction1Types_ConstantField_H

#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                        Class ConstantField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ConstantField
:
    public PatchFunction1<Type>
{
    // Private Data

        //- Is uniform?
        bool isUniform_;

        //- If uniform the uniformValue
        Type uniformValue_;

        //- ConstantField value
        Field<Type> value_;


    // Private Member Functions

        //- Helper to read value from dictionary entry
        static Field<Type> getValue
        (
            const word& entryName,
            const entry* eptr,  // primitive or dictionary entry
            const dictionary& dict,  // For error context or for content
            const label len,
            bool& isUniform,
            Type& uniformValue
        );

        //- No copy assignment
        void operator=(const ConstantField<Type>&) = delete;


public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a uniform value
        ConstantField
        (
            const polyPatch& pp,
            const word& entryName,
            const Type& uniformValue,
            const dictionary& dict = dictionary::null,
            const bool faceValues = true
        );

        //- Construct from components
        ConstantField
        (
            const polyPatch& pp,
            const word& entryName,
            const bool isUniform,
            const Type& uniformValue,
            const Field<Type>& fieldValues,
            const dictionary& dict = dictionary::null,
            const bool faceValues = true
        );

        //- Construct from entry name and dictionary
        ConstantField
        (
            const polyPatch& pp,
            const word& redirectType,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Construct from primitive entry, entry name and dictionary
        ConstantField
        (
            const polyPatch& pp,
            const entry* eptr,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct setting patch
        explicit ConstantField
        (
            const ConstantField<Type>& rhs,
            const polyPatch& pp
        );

        //- Copy construct
        explicit ConstantField(const ConstantField<Type>& rhs);

        //- Return a clone
        virtual tmp<PatchFunction1<Type>> clone() const
        {
            return PatchFunction1<Type>::Clone(*this);
        }

        //- Return a clone, setting the patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const
        {
            return PatchFunction1<Type>::Clone(*this, pp);
        }


    //- Destructor
    virtual ~ConstantField() = default;


    // Member Functions

        //- Value is independent of x
        virtual inline bool constant() const
        {
            return true;
        }

        //- Is value uniform (i.e. independent of coordinate)
        virtual inline bool uniform() const
        {
            return isUniform_ && PatchFunction1<Type>::uniform();
        }


        // Evaluation

            //- Return constant value
            virtual inline tmp<Field<Type>> value(const scalar x) const;

            //- Integrate between two values
            virtual inline tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<Type>& pf1,
                const labelList& addr
            );


        // I-O

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ConstantFieldI.H"

#ifdef NoRepository
    #include "ConstantField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::genericPatchFieldBase

Description
    Generic infrastructure for reading/writing unknown patch types.

SourceFiles
    genericPatchFieldBase.C
    genericPatchFieldBaseTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_genericPatchFieldBase_H
#define Foam_genericPatchFieldBase_H

#include "dictionary.H"
#include "primitiveFields.H"
#include "HashPtrTable.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class genericPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

class genericPatchFieldBase
{
    // Private Data

        //- The non-generic patch name
        word actualTypeName_;

        //- Copy of the input dictionary (used for transcribing)
        dictionary dict_;


    // The data fields read.
    // Do not need to handle generic for 'label' since these
    // will be limited to known field types.

        HashPtrTable<scalarField> scalarFields_;
        HashPtrTable<vectorField> vectorFields_;
        HashPtrTable<sphericalTensorField> sphTensorFields_;
        HashPtrTable<symmTensorField> symmTensorFields_;
        HashPtrTable<tensorField> tensorFields_;


    // Private Member Functions

        bool checkFieldSize
        (
            const label fieldSize,
            const label patchSize,
            const word& patchName,
            const keyType& key,
            const IOobject& io
        ) const;


protected:

    // Protected Member Functions

        //- Add error message to FatalError about solving with
        //- generic condition
        void genericFatalSolveError
        (
            const word& patchName,
            const IOobject& io
        ) const;

        //- FatalError for missing entry
        void reportMissingEntry
        (
            const word& entryName,
            const word& patchName,
            const IOobject& io
        ) const;

        void processGeneric
        (
            const label patchSize,
            const word& patchName,
            const IOobject& io,
            const bool separateValue
        );

        bool processEntry
        (
            const entry& dEntry,
            const label patchSize,
            const word& patchName,
            const IOobject& io
        );

        //- Write a single entry, with lookup of hashed values
        void putEntry(const entry& e, Ostream& os) const;

        //- Write all generic entries from dictionary,
        //- optionally treating the "value" entry separately
        void writeGeneric(Ostream& os, const bool separateValue) const;

        //- Implementation for construct with mapper
        template<class MapperType>
        void mapGeneric
        (
            const genericPatchFieldBase& rhs,
            const MapperType& mapper
        );

        //- Implementation for autoMap of self given a mapping object
        template<class MapperType>
        void autoMapGeneric(const MapperType& mapper);

        //- Implementation for reverse map given patch field onto this
        //- patch field
        void rmapGeneric
        (
            const genericPatchFieldBase& rhs,
            const labelList& addr
        );


    // Constructors

        //- Partial copy construct. Only copy type and dictionary
        genericPatchFieldBase(const Foam::zero, const genericPatchFieldBase&);


public:

    // Constructors

        //- Default construct, generally not useful.
        genericPatchFieldBase() = default;

        //- Copy construct
        genericPatchFieldBase(const genericPatchFieldBase&) = default;

        //- Move construct
        genericPatchFieldBase(genericPatchFieldBase&&) = default;

        //- Initialize from dictionary
        explicit genericPatchFieldBase(const dictionary& dict);


    // Member Functions

        //- The actual patch type
        const word& actualType() const noexcept
        {
            return actualTypeName_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "genericPatchFieldBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

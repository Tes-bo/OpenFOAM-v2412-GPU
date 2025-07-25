/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::DirectFieldMapper

Description
    A templated direct mapper for the given FieldMapper type

\*---------------------------------------------------------------------------*/

#ifndef Foam_directFieldMapper_H
#define Foam_directFieldMapper_H

#include "FieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class FieldMapperType> class DirectFieldMapper;

// Standard Types

//- A FieldMapper with direct mapping
typedef DirectFieldMapper<FieldMapper> directFieldMapper;


/*---------------------------------------------------------------------------*\
                      Class DirectFieldMapper Declaration
\*---------------------------------------------------------------------------*/

template<class FieldMapperType>
class DirectFieldMapper
:
    public FieldMapperType
{
    // Private Data

        //- Addressing from new back to old
        const labelUList& directAddressing_;

        //- Does map contain any unmapped values
        bool hasUnmapped_;


    // Private Member Functions

        //- Any negative (unmapped) values in the addressing?
        static bool hasUnmappedEntry(const labelUList& directAddr)
        {
            for (const label val : directAddr)
            {
                if (val < 0) return true;  // early exit
            }
            return false;
        }

public:

    // Public Types

        //- The base mapper type
        typedef FieldMapperType mapper_type;


    // Constructors

        //- Construct given addressing, check for unmapped (negative) values
        explicit DirectFieldMapper
        (
            const labelUList& directAddr,
            const bool checkUnmapped = true
        )
        :
            FieldMapperType(),
            directAddressing_(directAddr),
            hasUnmapped_(checkUnmapped && hasUnmappedEntry(directAddr))
        {}


    //- Destructor
    virtual ~DirectFieldMapper() = default;


    // Member Functions

        //- True if directAddressing is not the null object (unallocated)
        virtual bool hasDirectAddressing() const
        {
            return notNull(directAddressing_);
        }

        //- The mapper size is given by the size of the direct addressing
        virtual label size() const
        {
            return directAddressing_.size();
        }

        //- It is a direct mapper
        virtual bool direct() const
        {
            return true;
        }

        //- Any unmapped values?
        virtual bool hasUnmapped() const
        {
            return hasUnmapped_;
        }

        //- Allow modification
        virtual bool& hasUnmapped()
        {
            return hasUnmapped_;
        }

        //- Return the direct addressing values
        virtual const labelUList& directAddressing() const
        {
            return directAddressing_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

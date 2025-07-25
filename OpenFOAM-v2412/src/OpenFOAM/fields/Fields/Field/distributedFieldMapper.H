/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::DistributedFieldMapper

Description
    A templated field mapper
    with direct mapping from local or remote quantities.

\*---------------------------------------------------------------------------*/

#ifndef Foam_distributedFieldMapper_H
#define Foam_distributedFieldMapper_H

#include "directFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class DirectFieldMapperType> class DistributedFieldMapper;

// Standard Types

//- A directFieldMapper with distributed (with local or remote) quantities
typedef DistributedFieldMapper<directFieldMapper> distributedFieldMapper;


/*---------------------------------------------------------------------------*\
                   Class DistributedFieldMapper Declaration
\*---------------------------------------------------------------------------*/

template<class DirectFieldMapperType>
class DistributedFieldMapper
:
    public DirectFieldMapperType
{
    // Private Data

        //- The distributed map
        const mapDistributeBase& distMap_;


public:

    // Public Types

        //- The base direct mapper type
        typedef DirectFieldMapperType mapper_type;


    // Constructors

        //- Construct with addressing and distribute map
        DistributedFieldMapper
        (
            const labelUList& directAddr,
            const mapDistributeBase& distMap
        )
        :
            DirectFieldMapperType(directAddr),
            distMap_(distMap)
        {}


    //- Destructor
    virtual ~DistributedFieldMapper() = default;


    // Member Functions

        //- The mapper size is given by the direct addressing size (if valid)
        //- or from the distributeMap construct size.
        virtual label size() const
        {
            return
            (
                DirectFieldMapperType::hasDirectAddressing()
              ? DirectFieldMapperType::size()
              : distMap_.constructSize()
            );
        }

        //- It is a distributed mapper
        virtual bool distributed() const
        {
            return true;
        }

        //- Return the distribution map
        virtual const mapDistributeBase& distributeMap() const
        {
            return distMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

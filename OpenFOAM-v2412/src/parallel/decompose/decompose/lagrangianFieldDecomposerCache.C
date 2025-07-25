/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "lagrangianFieldDecomposer.H"

#include "labelIOField.H"
#include "labelFieldIOField.H"
#include "scalarIOField.H"
#include "scalarFieldIOField.H"
#include "vectorIOField.H"
#include "vectorFieldIOField.H"
#include "sphericalTensorIOField.H"
#include "sphericalTensorFieldIOField.H"
#include "symmTensorIOField.H"
#include "symmTensorFieldIOField.H"
#include "tensorIOField.H"
#include "tensorFieldIOField.H"

// * * * * * * * * * * * * * * * * Declarations  * * * * * * * * * * * * * * //

namespace Foam
{

// All lagrangian field/field-field types
class lagrangianFieldDecomposer::fieldsCache::privateCache
{
public:

    #undef  declareField
    #define declareField(Type)                                                \
    PtrList<PtrList<IOField<Type>>> Type##Fields_;                            \
    PtrList<PtrList<CompactIOField<Field<Type>, Type>>> Type##FieldFields_;

    declareField(label);
    declareField(scalar);
    declareField(vector);
    declareField(sphericalTensor);
    declareField(symmTensor);
    declareField(tensor);
    #undef declareField

    bool empty() const noexcept { return labelFields_.empty(); }

    label size() const noexcept { return labelFields_.size(); }

    void resize(const label len)
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            Type##Fields_.resize(len);                                \
            Type##FieldFields_.resize(len);                           \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }

    void readAll(const label cloudi, const IOobjectList& lagrangianObjects)
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            lagrangianFieldDecomposer::readFields                     \
            (                                                         \
                cloudi,                                               \
                lagrangianObjects,                                    \
                Type##Fields_                                         \
            );                                                        \
            lagrangianFieldDecomposer::readFieldFields                \
            (                                                         \
                cloudi,                                               \
                lagrangianObjects,                                    \
                Type##FieldFields_                                    \
            );                                                        \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }

    void decomposeAll
    (
        const label cloudi,
        const fileName& cloudDir,
        const lagrangianFieldDecomposer& decomposer,
        bool report  /* unused */
    ) const
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            decomposer.decomposeFields                                \
            (                                                         \
                cloudDir,                                             \
                Type##Fields_[cloudi]                                 \
            );                                                        \
            decomposer.decomposeFieldFields                           \
            (                                                         \
                cloudDir,                                             \
                Type##FieldFields_[cloudi]                            \
            );                                                        \
        }

        doLocalCode(label);
        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }
};

}  // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianFieldDecomposer::fieldsCache::fieldsCache()
:
    cache_(new privateCache)
{}


Foam::lagrangianFieldDecomposer::fieldsCache::fieldsCache
(
    const label nClouds
)
:
    cache_(new privateCache)
{
    cache_->resize(nClouds);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Destructor not in header (used incomplete type)
Foam::lagrangianFieldDecomposer::fieldsCache::~fieldsCache()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::lagrangianFieldDecomposer::fieldsCache::empty() const
{
    return (!cache_ || cache_->empty());
}


Foam::label Foam::lagrangianFieldDecomposer::fieldsCache::size() const
{
    return (cache_ ? cache_->size() : label(0));
}


void Foam::lagrangianFieldDecomposer::fieldsCache::clear()
{
    cache_.reset(new privateCache);
}


void Foam::lagrangianFieldDecomposer::fieldsCache::resize
(
    const label nClouds
)
{
    if (cache_)
    {
        cache_->resize(nClouds);
    }
}


void Foam::lagrangianFieldDecomposer::fieldsCache::readAllFields
(
    const label cloudi,
    const IOobjectList& lagrangianObjects
)
{
    if (cache_)
    {
        cache_->readAll(cloudi, lagrangianObjects);
    }
}


void Foam::lagrangianFieldDecomposer::fieldsCache::decomposeAllFields
(
    const label cloudi,
    const fileName& cloudDir,
    const lagrangianFieldDecomposer& decomposer,
    bool report
) const
{
    if (cache_)
    {
        cache_->decomposeAll(cloudi, cloudDir, decomposer, report);
    }
}


// ************************************************************************* //

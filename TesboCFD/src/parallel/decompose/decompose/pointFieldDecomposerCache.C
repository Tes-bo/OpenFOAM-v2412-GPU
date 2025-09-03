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

#include "pointFieldDecomposer.H"
#include "fieldsDistributor.H"
#include "pointFields.H"
#include "IOobjectList.H"
#include "PtrListOps.H"

// * * * * * * * * * * * * * * * * Declarations  * * * * * * * * * * * * * * //

namespace Foam
{

// All point field types
class pointFieldDecomposer::fieldsCache::privateCache
{
public:

    #undef  declareField
    #define declareField(Type)                                                \
    PtrList<GeometricField<Type, pointPatchField, pointMesh>> Type##Fields_;

    declareField(scalar);
    declareField(vector);
    declareField(sphericalTensor);
    declareField(symmTensor);
    declareField(tensor);
    #undef declareField

    label size() const noexcept
    {
        label count = 0;

        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            count += Type##Fields_.size();                            \
        }

        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);
        #undef doLocalCode

        return count;
    }

    bool empty() const noexcept { return !size(); }

    void readAll(const pointMesh& mesh, const IOobjectList& objects)
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            fieldsDistributor::readFields                             \
            (                                                         \
                mesh,                                                 \
                objects,                                              \
                Type##Fields_,                                        \
                false  /* readOldTime = false */                      \
            );                                                        \
        }

        doLocalCode(scalar);
        doLocalCode(vector);
        doLocalCode(sphericalTensor);
        doLocalCode(symmTensor);
        doLocalCode(tensor);

        #undef doLocalCode
    }

    template<class GeoField>
    static void decompose
    (
        const pointFieldDecomposer& decomposer,
        const PtrList<GeoField>& fields,
        bool report
    )
    {
        if (!fields.empty())
        {
            if (report)
            {
                Info<< "  "
                    << pTraits<typename GeoField::value_type>::typeName
                    << "s: "
                    << flatOutput(PtrListOps::names(fields)) << nl;
            }

            decomposer.decomposeFields(fields);
        }
    }

    void decomposeAll
    (
        const pointFieldDecomposer& decomposer,
        bool report
    ) const
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            decompose(decomposer, Type##Fields_, report);             \
        }

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

Foam::pointFieldDecomposer::fieldsCache::fieldsCache()
:
    cache_(new privateCache)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Destructor not in header (used incomplete type)
Foam::pointFieldDecomposer::fieldsCache::~fieldsCache()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pointFieldDecomposer::fieldsCache::empty() const
{
    return (!cache_ || cache_->empty());
}


Foam::label Foam::pointFieldDecomposer::fieldsCache::size() const
{
    return (cache_ ? cache_->size() : label(0));
}


void Foam::pointFieldDecomposer::fieldsCache::clear()
{
    cache_.reset(new privateCache);
}


void Foam::pointFieldDecomposer::fieldsCache::readAllFields
(
    const pointMesh& mesh,
    const IOobjectList& objects
)
{
    if (cache_)
    {
        cache_->readAll(mesh, objects);
    }
}


void Foam::pointFieldDecomposer::fieldsCache::decomposeAllFields
(
    const pointFieldDecomposer& decomposer,
    bool report
) const
{
    if (cache_)
    {
        cache_->decomposeAll(decomposer, report);
    }
}


// ************************************************************************* //

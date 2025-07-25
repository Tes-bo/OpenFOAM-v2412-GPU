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

#include "faFieldDecomposer.H"
#include "fieldsDistributor.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "IOobjectList.H"
#include "PtrListOps.H"

// * * * * * * * * * * * * * * * * Declarations  * * * * * * * * * * * * * * //

namespace Foam
{

// All area/edge field types
class faFieldDecomposer::fieldsCache::privateCache
{
public:

    #undef  declareField
    #define declareField(Type)                                                \
    PtrList<GeometricField<Type, faPatchField, areaMesh>> Type##AreaFields_;  \
    PtrList<GeometricField<Type, faePatchField, edgeMesh>> Type##EdgeFields_;

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
            count += Type##AreaFields_.size();                        \
            count += Type##EdgeFields_.size();                        \
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

    void readAll
    (
        const faMesh& mesh,
        const IOobjectList& objects
    )
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            fieldsDistributor::readFields                             \
            (                                                         \
                mesh,                                                 \
                objects,                                              \
                Type##AreaFields_                                     \
            );                                                        \
            fieldsDistributor::readFields                             \
            (                                                         \
                mesh,                                                 \
                objects,                                              \
                Type##AreaFields_,                                    \
                false  /* readOldTime = false */                      \
            );                                                        \
            fieldsDistributor::readFields                             \
            (                                                         \
                mesh,                                                 \
                objects,                                              \
                Type##EdgeFields_,                                    \
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

    template<class BoolListType>
    void readAll
    (
        const BoolListType& haveMeshOnProc,
        const faMeshSubset*& subsetter,
        const faMesh& mesh,
        IOobjectList& objects
    )
    {
        #undef  doLocalCode
        #define doLocalCode(Type)                                     \
        {                                                             \
            fieldsDistributor::readFields                             \
            (                                                         \
                haveMeshOnProc,                                       \
                subsetter,                                            \
                mesh,                                                 \
                objects,                                              \
                Type##AreaFields_                                     \
            );                                                        \
            fieldsDistributor::readFields                             \
            (                                                         \
                haveMeshOnProc,                                       \
                subsetter,                                            \
                mesh,                                                 \
                objects,                                              \
                Type##EdgeFields_                                     \
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
        const faFieldDecomposer& decomposer,
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

    void decomposeAll(const faFieldDecomposer& decomposer, bool report) const
    {
        #undef  doLocalCode
        #define doLocalCode(Flavour)                                          \
        {                                                                     \
            decompose(decomposer, scalar##Flavour##Fields_, report);          \
            decompose(decomposer, vector##Flavour##Fields_, report);          \
            decompose(decomposer, sphericalTensor##Flavour##Fields_, report); \
            decompose(decomposer, symmTensor##Flavour##Fields_, report);      \
            decompose(decomposer, tensor##Flavour##Fields_, report);          \
        }

        doLocalCode(Area);
        doLocalCode(Edge);

        #undef doLocalCode
    }
};

}  // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faFieldDecomposer::fieldsCache::fieldsCache()
:
    cache_(new privateCache)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Destructor not in header (used incomplete type)
Foam::faFieldDecomposer::fieldsCache::~fieldsCache()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faFieldDecomposer::fieldsCache::empty() const
{
    return (!cache_ || cache_->empty());
}


Foam::label Foam::faFieldDecomposer::fieldsCache::size() const
{
    return (cache_ ? cache_->size() : label(0));
}


void Foam::faFieldDecomposer::fieldsCache::clear()
{
    cache_.reset(new privateCache);
}


void Foam::faFieldDecomposer::fieldsCache::readAllFields
(
    const faMesh& mesh,
    const IOobjectList& objects
)
{
    if (cache_)
    {
        cache_->readAll(mesh, objects);
    }
}


void Foam::faFieldDecomposer::fieldsCache::readAllFields
(
    const bitSet& haveMeshOnProc,
    const faMeshSubset* subsetter,
    const faMesh& mesh,
    IOobjectList& objects
)
{
    if (cache_)
    {
        cache_->readAll(haveMeshOnProc, subsetter, mesh, objects);
    }
}


void Foam::faFieldDecomposer::fieldsCache::readAllFields
(
    const boolList& haveMeshOnProc,
    const faMeshSubset* subsetter,
    const faMesh& mesh,
    IOobjectList& objects
)
{
    if (cache_)
    {
        cache_->readAll(haveMeshOnProc, subsetter, mesh, objects);
    }
}


void Foam::faFieldDecomposer::fieldsCache::decomposeAllFields
(
    const faFieldDecomposer& decomposer,
    bool report
) const
{
    if (cache_)
    {
        cache_->decomposeAll(decomposer, report);
    }
}


// ************************************************************************* //

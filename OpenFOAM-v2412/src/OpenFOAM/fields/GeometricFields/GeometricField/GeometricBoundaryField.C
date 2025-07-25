/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "GeometricBoundaryField.H"
#include "globalMeshData.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
template<class CheckPatchFieldType>
bool Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::checkConsistency
(
    const scalar tol,
    const bool doExit
) const
{
    if (!this->size())
    {
        return true;
    }

    if (debug&2)
    {
        const auto& pfld0 = this->operator[](0);
        PoutInFunction
            << " Checking boundary consistency for field "
            << pfld0.internalField().name()
            << endl;
    }

    auto& bfld = const_cast<GeometricBoundaryField<Type, PatchField, GeoMesh>&>
    (
        *this
    );


    // Store old value
    List<Field<Type>> oldBfld(this->size());
    boolList oldUpdated(this->size());
    //Note: areaFields (finiteArea) do not have manipulatedMatrix() flag. TBD.
    //boolList oldManipulated(this->size());

    for (auto& pfld : bfld)
    {
        if (isA<CheckPatchFieldType>(pfld))
        {
            const label patchi = pfld.patch().index();
            oldUpdated[patchi] = pfld.updated();
            oldBfld[patchi] = pfld;
            //oldManipulated[patchi] = pfld.manipulatedMatrix();
        }
    }


    // Re-evaluate
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : bfld)
        {
            if (isA<CheckPatchFieldType>(pfld))
            {
                pfld.initEvaluate(UPstream::commsTypes::nonBlocking);
            }
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : bfld)
        {
            if (isA<CheckPatchFieldType>(pfld))
            {
                pfld.evaluate(UPstream::commsTypes::nonBlocking);
            }
        }
    }


    // Check
    bool ok = true;
    for (auto& pfld : bfld)
    {
        if (isA<CheckPatchFieldType>(pfld))
        {
            const label patchi = pfld.patch().index();
            const auto& oldPfld = oldBfld[patchi];

            forAll(pfld, facei)
            {
                if (mag(pfld[facei]-oldPfld[facei]) > tol)
                {
                    ok = false;
                    break;
                }
            }

            if (!ok)
            {
                if (doExit)
                {
                    FatalErrorInFunction << "Field "
                        << pfld.internalField().name()
                        << " is not evaluated?"
                        << " On patch " << pfld.patch().name()
                        << " type " << pfld.type()
                        << " : average of field = "
                        << average(oldPfld)
                        << ". Average of evaluated field = "
                        << average(pfld)
                        << ". Difference:" << average(pfld-oldPfld)
                        << ". Tolerance:" << tol
                        << exit(FatalError);
                }
                else
                {
                    WarningInFunction << "Field "
                        << pfld.internalField().name()
                        << " is not evaluated?"
                        << " On patch " << pfld.patch().name()
                        << " type " << pfld.type()
                        << " : average of field = "
                        << average(oldPfld)
                        << ". Average of evaluated field = "
                        << average(pfld)
                        << ". Difference:" << average(pfld-oldPfld)
                        << ". Tolerance:" << tol
                        << endl;

                    // Skip other patches
                    break;
                }
            }
        }
    }

    // Restore bfld, updated
    for (auto& pfld : bfld)
    {
        if (isA<CheckPatchFieldType>(pfld))
        {
            const label patchi = pfld.patch().index();
            pfld.setUpdated(oldUpdated[patchi]);
            Field<Type>& vals = pfld;
            vals = std::move(oldBfld[patchi]);
            //pfld.setManipulated(oldManipulated[patchi]);
        }
    }

    if (debug&2)
    {
        const auto& pfld0 = this->operator[](0);
        PoutInFunction
            << " Result of checking for field "
            << pfld0.internalField().name() << " : " << ok << endl;
    }

    return ok;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::readField
(
    const Internal& iField,
    const dictionary& dict
)
{
    // DebugInFunction << nl;

    // Clear the boundary field if already initialised
    this->clear();

    this->resize(bmesh_.size());

    label nUnset = this->size();

    // 1. Handle explicit patch names. Note that there can be only one explicit
    //    patch name since is key of dictionary.

    for (const entry& dEntry : dict)
    {
        const auto* subdict = dEntry.dictPtr();

        if (subdict && dEntry.keyword().isLiteral())
        {
            const label patchi = bmesh_.findPatchID(dEntry.keyword());

            if (patchi != -1)
            {
                this->set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        bmesh_[patchi],
                        iField,
                        *subdict
                    )
                );
                --nUnset;
            }
        }
    }

    if (nUnset == 0)
    {
        return;
    }


    // 2. Patch-groups. (using non-wild card entries of dictionaries)
    // (patchnames already matched above)
    // Note: in reverse order of entries in the dictionary (last
    // patchGroups wins). This is so it is consistent with dictionary wildcard
    // behaviour
    for (auto iter = dict.crbegin(); iter != dict.crend(); ++iter)
    {
        const entry& dEntry = *iter;
        const auto* subdict = dEntry.dictPtr();

        if (subdict && dEntry.keyword().isLiteral())
        {
            const labelList patchIds =
                bmesh_.indices(dEntry.keyword(), true); // use patchGroups

            for (const label patchi : patchIds)
            {
                if (!this->set(patchi))
                {
                    this->set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            bmesh_[patchi],
                            iField,
                            *subdict
                        )
                    );
                }
            }
        }
    }


    // 3. Wildcard patch overrides
    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            if (bmesh_[patchi].type() == emptyPolyPatch::typeName)
            {
                this->set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        emptyPolyPatch::typeName,
                        bmesh_[patchi],
                        iField
                    )
                );
            }
            else
            {
                const auto* subdict = dict.findDict(bmesh_[patchi].name());

                if (subdict)
                {
                    this->set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            bmesh_[patchi],
                            iField,
                            *subdict
                        )
                    );
                }
            }
        }
    }


    // Check for any unset patches
    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            if (bmesh_[patchi].type() == cyclicPolyPatch::typeName)
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find patchField entry for cyclic "
                    << bmesh_[patchi].name() << endl
                    << "Is your field uptodate with split cyclics?" << endl
                    << "Run foamUpgradeCyclics to convert mesh and fields"
                    << " to split cyclics." << exit(FatalIOError);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find patchField entry for "
                    << bmesh_[patchi].name() << exit(FatalIOError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    // DebugInFunction << nl;

    forAll(bmesh_, patchi)
    {
        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchi],
                iField
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const wordList& patchFieldTypes,
    const wordList& constraintTypes
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    // DebugInFunction << nl;

    if
    (
        patchFieldTypes.size() != this->size()
     || (constraintTypes.size() && (constraintTypes.size() != this->size()))
    )
    {
        FatalErrorInFunction
            << "Incorrect number of patch type specifications given" << nl
            << "    Number of patches in mesh = " << bmesh.size()
            << " number of patch type specifications = "
            << patchFieldTypes.size()
            << abort(FatalError);
    }

    if (constraintTypes.size())
    {
        forAll(bmesh_, patchi)
        {
            this->set
            (
                patchi,
                PatchField<Type>::New
                (
                    patchFieldTypes[patchi],
                    constraintTypes[patchi],
                    bmesh_[patchi],
                    iField
                )
            );
        }
    }
    else
    {
        forAll(bmesh_, patchi)
        {
            this->set
            (
                patchi,
                PatchField<Type>::New
                (
                    patchFieldTypes[patchi],
                    bmesh_[patchi],
                    iField
                )
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const PtrList<PatchField<Type>>& ptfl
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    // DebugInFunction << nl;

    forAll(bmesh_, patchi)
    {
        this->set(patchi, ptfl[patchi].clone(iField));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const Internal& iField,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    // DebugInFunction << nl;

    forAll(bmesh_, patchi)
    {
        this->set(patchi, btf[patchi].clone(iField));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const Internal& iField,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf,
    const labelList& patchIDs,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    // DebugInFunction << nl;

    for (const label patchi : patchIDs)
    {
        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchi],
                iField
            )
        );
    }

    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            this->set(patchi, btf[patchi].clone(iField));
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf
)
:
    FieldField<PatchField, Type>(btf),
    bmesh_(btf.bmesh_)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const Internal& iField,
    const dictionary& dict
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    readField(iField, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::updateCoeffs()
{
    // DebugInFunction << nl;

    for (auto& pfld : *this)
    {
        pfld.updateCoeffs();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluate
(
    const UPstream::commsTypes commsType
)
{
    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            pfld.initEvaluate(commsType);
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : *this)
        {
            pfld.evaluate(commsType);
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (schedEval.init)
            {
                pfld.initEvaluate(commsType);
            }
            else
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class UnaryPredicate>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluate_if
(
    const UnaryPredicate& pred,
    const UPstream::commsTypes commsType
)
{
    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            if (pred(pfld))
            {
                pfld.initEvaluate(commsType);
            }
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : *this)
        {
            if (pred(pfld))
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (pred(pfld))
            {
                if (schedEval.init)
                {
                    pfld.initEvaluate(commsType);
                }
                else
                {
                    pfld.evaluate(commsType);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateSelected
(
    const labelUList& patchIDs
)
{
    const label startOfRequests = UPstream::nRequests();

    for (const label patchi : patchIDs)
    {
        auto& pfld = (*this)[patchi];

        DebugInfo<< "Updating " << pfld.patch().name() << endl;

        pfld.initEvaluate(UPstream::commsTypes::nonBlocking);
    }

    // Wait for outstanding requests (non-blocking)
    UPstream::waitRequests(startOfRequests);

    for (const label patchi : patchIDs)
    {
        auto& pfld = (*this)[patchi];

        pfld.evaluate(UPstream::commsTypes::nonBlocking);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateLocal
(
    const UPstream::commsTypes commsType
)
{
    // DebugInFunction << nl;

    if (!localConsistency)
    {
        return;
    }

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            pfld.initEvaluateLocal(commsType);
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : *this)
        {
            pfld.evaluateLocal(commsType);
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (schedEval.init)
            {
                pfld.initEvaluateLocal(commsType);
            }
            else
            {
                pfld.evaluateLocal(commsType);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class CoupledPatchType>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateCoupled
(
    const UPstream::commsTypes commsType
)
{
    // Alternative (C++14)
    //
    // this->evaluate_if
    // (
    //     [](const auto& pfld) -> bool
    //     {
    //         const auto* cpp = isA<CoupledPatchType>(pfld.patch());
    //         return (cpp && cpp->coupled());
    //     },
    //     commsType
    // );

    // DebugInFunction << nl;

    if
    (
        commsType == UPstream::commsTypes::buffered
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            const auto* cpp = isA<CoupledPatchType>(pfld.patch());

            if (cpp && cpp->coupled())
            {
                pfld.initEvaluate(commsType);
            }
        }

        // Wait for outstanding requests (non-blocking)
        UPstream::waitRequests(startOfRequests);

        for (auto& pfld : *this)
        {
            const auto* cpp = isA<CoupledPatchType>(pfld.patch());

            if (cpp && cpp->coupled())
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            const auto* cpp = isA<CoupledPatchType>(pfld.patch());

            if (cpp && cpp->coupled())
            {
                if (schedEval.init)
                {
                    pfld.initEvaluate(commsType);
                }
                else
                {
                    pfld.evaluate(commsType);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type " << int(commsType) << nl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::types() const
{
    const FieldField<PatchField, Type>& pff = *this;

    wordList list(pff.size());

    forAll(pff, patchi)
    {
        list[patchi] = pff[patchi].type();
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::
boundaryInternalField() const
{
    auto tresult = tmp<GeometricBoundaryField<Type, PatchField, GeoMesh>>::New
    (
        DimensionedField<Type, GeoMesh>::null(),
        *this
    );

    auto& result = tresult.ref();

    forAll(result, patchi)
    {
        result[patchi] == this->operator[](patchi).patchInternalField();
    }

    return tresult;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::LduInterfaceFieldPtrsList<Type>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::interfaces() const
{
    LduInterfaceFieldPtrsList<Type> list(this->size());

    forAll(list, patchi)
    {
        list.set
        (
            patchi,
            isA<LduInterfaceField<Type>>(this->operator[](patchi))
        );
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::lduInterfaceFieldPtrsList
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::
scalarInterfaces() const
{
    lduInterfaceFieldPtrsList list(this->size());

    forAll(list, patchi)
    {
        list.set
        (
            patchi,
            isA<lduInterfaceField>(this->operator[](patchi))
        );
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);
    this->writeEntries(os);
    os.endBlock();

    os.check(FUNCTION_NAME);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::writeEntries
(
    Ostream& os
) const
{
    for (const auto& pfld : *this)
    {
        os.beginBlock(pfld.patch().name());
        os << pfld;
        os.endBlock();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
bool Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::check
(
) const
{
    // Dummy op - template specialisations provide logic (usually call
    // to checkConsistency)
    return true;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const FieldField<PatchField, Type>& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const Type& val
)
{
    FieldField<PatchField, Type>::operator=(val);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const FieldField<PatchField, Type>& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const Type& val
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == val;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    os << static_cast<const FieldField<PatchField, Type>&>(bf);
    return os;
}


// ************************************************************************* //

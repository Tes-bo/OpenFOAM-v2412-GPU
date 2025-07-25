/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "SemiImplicitSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmSup.H"
#include "Constant.H"
#include "Tuple2.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Type>
const Foam::Enum
<
    typename Foam::fv::SemiImplicitSource<Type>::volumeModeType
>
Foam::fv::SemiImplicitSource<Type>::volumeModeTypeNames_
({
    { volumeModeType::vmAbsolute, "absolute" },
    { volumeModeType::vmSpecific, "specific" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::setFieldCoeffs
(
    const dictionary& dict,
    const word& keyExplicit,
    const word& keyImplicit
)
{
    label count = dict.size();

    fieldNames_.resize_nocopy(count);

    Su_.clear();
    Sp_.clear();
    Su_.reserve(count);
    Sp_.reserve(count);

    driverSu_.clear();
    driverSp_.clear();
    driverSu_.reserve(count);
    driverSp_.reserve(count);

    valueExprSu_.clear();
    valueExprSp_.clear();
    valueExprSu_.reserve(count);
    valueExprSp_.reserve(count);

    fv::option::resetApplied();

    word modelType;
    Tuple2<Type, scalar> sourceRates;

    count = 0;

    for (const entry& dEntry : dict)
    {
        const word& fieldName = dEntry.keyword();
        bool ok = false;

        if (dEntry.isDict())
        {
            const dictionary& subDict = dEntry.dict();

            const entry* eptr;

            if
            (
                (eptr = subDict.findEntry(keyExplicit, keyType::LITERAL))
             != nullptr
            )
            {
                ok = true;

                if
                (
                    eptr->isDict()
                 && eptr->dict().readEntry("type", modelType, keyType::LITERAL)
                 && (modelType == "exprField")
                )
                {
                    const dictionary& exprDict = eptr->dict();

                    valueExprSu_.emplace_set(fieldName);
                    valueExprSu_[fieldName].readEntry("expression", exprDict);

                    driverSu_.set
                    (
                        fieldName,
                        new expressions::volumeExprDriver(mesh_, exprDict)
                    );
                }
                else
                {
                    Su_.set
                    (
                        fieldName,
                        Function1<Type>::New(keyExplicit, subDict, &mesh_)
                    );
                }
            }

            if
            (
                (eptr = subDict.findEntry(keyImplicit, keyType::LITERAL))
             != nullptr
            )
            {
                ok = true;

                if
                (
                    eptr->isDict()
                 && eptr->dict().readEntry("type", modelType, keyType::LITERAL)
                 && (modelType == "exprField")
                )
                {
                    const dictionary& exprDict = eptr->dict();

                    valueExprSp_.emplace_set(fieldName);
                    valueExprSp_[fieldName].readEntry("expression", exprDict);

                    driverSp_.set
                    (
                        fieldName,
                        new expressions::volumeExprDriver(mesh_, exprDict)
                    );
                }
                else
                {
                    Sp_.set
                    (
                        fieldName,
                        Function1<scalar>::New(keyImplicit, subDict, &mesh_)
                    );
                }
            }
        }
        else
        {
            // Non-dictionary form

            dEntry.readEntry(sourceRates);

            ok = true;

            Su_.set
            (
                fieldName,
                new Function1Types::Constant<Type>
                (
                    keyExplicit,
                    sourceRates.first()
                )
            );
            Sp_.set
            (
                fieldName,
                new Function1Types::Constant<scalar>
                (
                    keyImplicit,
                    sourceRates.second()
                )
            );
        }

        if (!ok)
        {
            FatalIOErrorInFunction(dict)
                << "Require at least one of "
                << keyExplicit << '/' << keyImplicit << " entries for "
                << "field: " << fieldName << endl
                << exit(FatalIOError);
        }

        fieldNames_[count] = fieldName;
        ++count;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fv::SemiImplicitSource<Type>::SemiImplicitSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    volumeMode_(vmAbsolute),
    VDash_(1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    return this->addSup(volScalarField::null(), eqn, fieldi);
}


template<class Type>
void Foam::fv::SemiImplicitSource<Type>::addSup
(
    const volScalarField& rho,
    fvMatrix<Type>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< "SemiImplicitSource<" << pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    const GeometricField<Type, fvPatchField, volMesh>& psi = eqn.psi();

    // Note: field name may deviate from psi name
    const word& fieldName = fieldNames_[fieldi];

    const scalar tmVal = mesh_.time().timeOutputValue();

    const dimensionSet SuDims(eqn.dimensions()/dimVolume);
    const dimensionSet SpDims(SuDims/psi.dimensions());

    // Explicit source
    {
        const auto iter1 = valueExprSu_.cfind(fieldName);
        const auto iter2 = Su_.cfind(fieldName);

        tmp<DimensionedField<Type, volMesh>> tsu;

        if (iter1.good())
        {
            const auto& valueExpr = iter1.val();

            typedef
                GeometricField<Type, fvPatchField, volMesh>
                ExprResultType;

            if (debug)
            {
                Info<< "Explicit expression source:" << nl
                    << ">>>>" << nl
                    << valueExpr.c_str() << nl
                    << "<<<<" << nl;
            }

            auto& driver = *(driverSu_[fieldName]);

            driver.clearVariables();

            if (notNull(rho))
            {
                driver.addContextObject("rho", &rho);
            }

            // Switch dimension checking off
            const bool oldDimChecking = dimensionSet::checking(false);

            driver.parse(valueExpr);

            // Restore dimension checking
            dimensionSet::checking(oldDimChecking);

            const ExprResultType* ptr =
                driver.template isResultType<ExprResultType>();

            if (!ptr)
            {
                FatalErrorInFunction
                    << "Expression for Su " << fieldName
                    << " evaluated to <" << driver.resultType()
                    << "> but expected <" << ExprResultType::typeName
                    << ">" << endl
                    << exit(FatalError);
            }
            else if (ptr->size() != mesh_.nCells())
            {
                FatalErrorInFunction
                    << "Expression for Su " << fieldName
                    << " evaluated to " << ptr->size()
                    << " instead of " << mesh_.nCells() << " values" << endl
                    << exit(FatalError);
            }

            if (notNull(rho))
            {
                driver.removeContextObject(&rho);
            }

            const Field<Type>& exprFld = ptr->primitiveField();

            tsu = DimensionedField<Type, volMesh>::New
            (
                name_ + fieldName + "Su",
                IOobject::NO_REGISTER,
                mesh_,
                dimensioned<Type>(SuDims, Zero)
            );

            if (this->useSubMesh())
            {
                for (const label celli : cells_)
                {
                    tsu.ref()[celli] = exprFld[celli]/VDash_;
                }
            }
            else
            {
                tsu.ref().field() = exprFld;

                if (!equal(VDash_, 1))
                {
                    tsu.ref().field() /= VDash_;
                }
            }
        }
        else if (iter2.good() && iter2.val()->good())
        {
            const dimensioned<Type> SuValue
            (
                "Su",
                SuDims,
                iter2.val()->value(tmVal)/VDash_
            );

            if (mag(SuValue.value()) <= ROOTVSMALL)
            {
                // No-op
            }
            else if (this->useSubMesh())
            {
                tsu = DimensionedField<Type, volMesh>::New
                (
                    name_ + fieldName + "Su",
                    IOobject::NO_REGISTER,
                    mesh_,
                    dimensioned<Type>(SuDims, Zero)
                );
                UIndirectList<Type>(tsu.ref(), cells_) = SuValue.value();
            }
            else
            {
                eqn += SuValue;
            }
        }

        if (tsu.valid())
        {
            eqn += tsu;
        }
    }


    // Implicit source
    {
        const auto iter1 = valueExprSp_.cfind(fieldName);
        const auto iter2 = Sp_.cfind(fieldName);

        tmp<DimensionedField<scalar, volMesh>> tsp;

        if (iter1.good())
        {
            const auto& valueExpr = iter1.val();

            typedef volScalarField ExprResultType;

            if (debug)
            {
                Info<< "Implicit expression source:" << nl
                    << ">>>>" << nl
                    << valueExpr.c_str() << nl
                    << "<<<<" << nl;
            }

            auto& driver = *(driverSp_[fieldName]);

            driver.clearVariables();

            if (notNull(rho))
            {
                driver.addContextObject("rho", &rho);
            }

            // Switch dimension checking off
            const bool oldDimChecking = dimensionSet::checking(false);

            driver.parse(valueExpr);

            // Restore dimension checking
            dimensionSet::checking(oldDimChecking);

            const ExprResultType* ptr =
                driver.template isResultType<ExprResultType>();

            if (!ptr)
            {
                FatalErrorInFunction
                    << "Expression for Sp " << fieldName
                    << " evaluated to <" << driver.resultType()
                    << "> but expected <" << ExprResultType::typeName
                    << ">" << endl
                    << exit(FatalError);
            }
            else if (ptr->size() != mesh_.nCells())
            {
                FatalErrorInFunction
                    << "Expression for Sp " << fieldName
                    << " evaluated to " << ptr->size()
                    << " instead of " << mesh_.nCells() << " values" << endl
                    << exit(FatalError);
            }

            if (notNull(rho))
            {
                driver.removeContextObject(&rho);
            }

            const Field<scalar>& exprFld = ptr->primitiveField();

            tsp = DimensionedField<scalar, volMesh>::New
            (
                name_ + fieldName + "Sp",
                IOobject::NO_REGISTER,
                mesh_,
                dimensioned<scalar>(SpDims, Zero)
            );

            if (this->useSubMesh())
            {
                for (const label celli : cells_)
                {
                    tsp.ref()[celli] = exprFld[celli]/VDash_;
                }
            }
            else
            {
                tsp.ref().field() = exprFld;

                if (!equal(VDash_, 1))
                {
                    tsp.ref().field() /= VDash_;
                }
            }
        }
        else if (iter2.good() && iter2.val()->good())
        {
            const dimensioned<scalar> SpValue
            (
                "Sp",
                SpDims,
                iter2.val()->value(tmVal)/VDash_
            );

            if (mag(SpValue.value()) <= ROOTVSMALL)
            {
                // No-op
            }
            else if (this->useSubMesh())
            {
                tsp = DimensionedField<scalar, volMesh>::New
                (
                    name_ + fieldName + "Sp",
                    IOobject::NO_REGISTER,
                    mesh_,
                    dimensioned<scalar>(SpDims, Zero)
                );
                UIndirectList<scalar>(tsp.ref(), cells_) = SpValue.value();
            }
            else
            {
                eqn += fvm::SuSp(SpValue, psi);
            }
        }

        if (tsp.valid())
        {
            eqn += fvm::SuSp(tsp, psi);
        }
    }
}


template<class Type>
bool Foam::fv::SemiImplicitSource<Type>::read(const dictionary& dict)
{
    VDash_ = 1;

    if (fv::cellSetOption::read(dict))
    {
        volumeMode_ = volumeModeTypeNames_.get("volumeMode", coeffs_);

        // Set volume normalisation
        if (volumeMode_ == vmAbsolute)
        {
            VDash_ = V_;
        }

        // Compatibility (2112 and earlier)
        const dictionary* injectDict =
            coeffs_.findDict("injectionRateSuSp", keyType::LITERAL);

        if (injectDict)
        {
            setFieldCoeffs
            (
                *injectDict,
                "Su",  // Su = explicit
                "Sp"   // Sp = implicit
            );
        }
        else
        {
            setFieldCoeffs
            (
                coeffs_.subDict("sources", keyType::LITERAL),
                "explicit",
                "implicit"
            );
        }

        return true;
    }

    return false;
}


// ************************************************************************* //

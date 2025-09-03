/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::expressions::exprDriver::exprDriver::weightedAverage
(
    const scalarField& wfield,
    const Field<Type>& fld
)
{
    if (isNull(wfield))
    {
        const label n = returnReduce(fld.size(), sumOp<label>());

        // stabilize
        if (!n)
        {
            return Zero;
        }

        return gSum(fld) / scalar(n);
    }

    // #ifdef FULLDEBUG
    // checkSize(wfield, fld);
    // #endif

    const scalar s = gSum(wfield);

    // stabilize
    if (mag(s) < ROOTVSMALL)
    {
        return Zero;
    }

    return gSum(wfield*fld) / s;
}


template<class Type>
Type Foam::expressions::exprDriver::exprDriver::weightedSum
(
    const scalarField& wfield,
    const Field<Type>& fld
)
{
    if (isNull(wfield))
    {
        return gSum(fld);
    }

    // #ifdef FULLDEBUG
    // checkSize(wfield, fld);
    // #endif

    return gSum(wfield*fld);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprDriver::getResult(bool wantPointData)
{
    if (!result_.isPointData(wantPointData))
    {
        FatalErrorInFunction
            << "Expected a" << (wantPointData ? " point" : "")
            << " field,  but found a" << (!wantPointData ? " point" : "")
            << " field" << nl
            << exit(FatalError);
    }

    return result_.getResult<Type>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::Function1<Type>* Foam::expressions::exprDriver::getFunction1Ptr
(
    const word& name,
    const HashTable<refPtr<Function1<Type>>>& tbl,
    wordList* listFailure
)
{
    const Function1<Type>* func = nullptr;

    const auto iter = tbl.cfind(name);

    if (iter.good())
    {
        func = iter.val().get();
    }

    if (!func && listFailure)
    {
        *listFailure = tbl.sortedToc();
    }

    return func;
}


template<class Type>
bool Foam::expressions::exprDriver::isFunction(const word& name) const
{
    // Currently only scalar, vector
    #undef doLocalCode
    #define doLocalCode(WhichType, MapperMember)                    \
    if (std::is_same<Type, WhichType>::value)                       \
    {                                                               \
        return bool                                                 \
        (                                                           \
            this->template getFunction1Ptr<WhichType>               \
            (                                                       \
                name, MapperMember                                  \
            )                                                       \
        );                                                          \
    }

    doLocalCode(scalar, scalarFuncs_);
    doLocalCode(vector, vectorFuncs_);
    #undef doLocalCode

    return false;
}


template<class Type>
Type Foam::expressions::exprDriver::getFunctionValue
(
    const word& name,
    const scalar x
) const
{
    const Function1<Type>* func = nullptr;

    wordList failed;

    do
    {
        // Currently only scalar, vector
        #undef doLocalCode
        #define doLocalCode(WhichType, MapperMember)                \
        if (std::is_same<Type, WhichType>::value)                   \
        {                                                           \
            const Function1<WhichType>* ptr =                       \
                this->template getFunction1Ptr<WhichType>           \
                (                                                   \
                    name, MapperMember, &failed                     \
                );                                                  \
            func = reinterpret_cast<const Function1<Type>*>(ptr);   \
            break;                                                  \
        }

        doLocalCode(scalar, scalarFuncs_);
        doLocalCode(vector, vectorFuncs_);
        #undef doLocalCode
    }
    while (false);

    // Error handling
    if (!failed.empty())
    {
        FatalErrorInFunction
            << "No mapping '" << name << " (" << pTraits<Type>::typeName
            << ") found." << nl
            << "Valid entries: "
            << flatOutput(failed) << nl
            << exit(FatalError);
    }

    if (func)
    {
        return func->value(x);
    }

    return pTraits<Type>::zero;
}


template<class Type>
void Foam::expressions::exprDriver::fillFunctionValues
(
    Field<Type>& result,
    const word& name,
    const scalarField& input
) const
{
    // #ifdef FULLDEBUG
    // checkSize(result, input);
    // #endif

    const Function1<Type>* func = nullptr;

    wordList failed;

    do
    {
        // Currently only scalar, vector
        #undef doLocalCode
        #define doLocalCode(WhichType, MapperMember)                \
        if (std::is_same<Type, WhichType>::value)                   \
        {                                                           \
            const Function1<WhichType>* ptr =                       \
                this->template getFunction1Ptr<WhichType>           \
                (                                                   \
                    name, MapperMember, &failed                     \
                );                                                  \
            func = reinterpret_cast<const Function1<Type>*>(ptr);   \
            break;                                                  \
        }

        doLocalCode(scalar, scalarFuncs_);
        doLocalCode(vector, vectorFuncs_);
        #undef doLocalCode
    }
    while (false);

    // Error handling
    if (!failed.empty())
    {
        FatalErrorInFunction
            << "No mapping '" << name << " (" << pTraits<Type>::typeName
            << ") found." << nl
            << "Valid entries: "
            << flatOutput(failed) << nl
            << exit(FatalError);
    }

    if (func)
    {
        const label len = min(result.size(), input.size());

        for (label i = 0; i < len; ++i)
        {
            result[i] = func->value(input[i]);
        }

        // Safety
        for (label i = len; i < result.size(); ++i)
        {
            result[i] = Zero;
        }

        return;
    }

    result = Zero;
}


template<class Type>
bool Foam::expressions::exprDriver::isLocalVariable
(
    const word& name,
    bool wantPointData,
    label expectedSize
) const
{
    DebugInfo
        << "Looking for local" << (wantPointData ? " point" : "")
        << " field name:" << name << " type:"
        << pTraits<Type>::typeName << " size:" << expectedSize;


    bool good = hasVariable(name);

    if (good)
    {
        const exprResult& var = variable(name);

        DebugInfo
            << " - found (" << var.valueType()
            << (var.isPointData() ? " point" : "") << ')';

        good = (var.isType<Type>() && var.isPointData(wantPointData));

        // Do size checking if requested
        if (good && expectedSize >= 0)
        {
            good = returnReduceAnd(var.size() == expectedSize);

            if (debug && !good)
            {
                Info<< " size is";
            }
        }
    }

    DebugInfo << (good ? " good" : " bad") << endl;

    return good;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprDriver::newField
(
    const Type& val
) const
{
    return tmp<Field<Type>>::New(size(), val);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::expressions::exprDriver::newPointField
(
    const Type& val
) const
{
    return tmp<Field<Type>>::New(pointSize(), val);
}


// ************************************************************************* //

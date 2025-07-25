/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "FunctionObjectValue.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1Types::FunctionObjectValue<Type>::value
(
    const scalar /*unused*/
) const
{
    const auto& props = this->time().functionObjects().propsDict();

    if (!props.hasResultObject(foName_))
    {
        if (haveDefaultValue_)
        {
            DebugInfo
                << "    Function object " << foName_
                << " not found; using default value " << defaultValue_
                << endl;

            return defaultValue_;
        }

        FatalErrorInFunction
            << "Function object " << foName_
            << " results not found. Valid objects with results include: "
            << props.objectResultNames()
            << exit(FatalError);
    }

    if (!props.hasResultObjectEntry(foName_, foResultName_))
    {
        if (haveDefaultValue_)
        {
            DebugInfo
                << "    Function object " << foName_
                << " result " << foResultName_
                << " not found; using default value " << defaultValue_
                << endl;

            return defaultValue_;
        }

        FatalErrorInFunction
            << "Function object " << foName_
            << " does not have a result field " << foResultName_ << nl
            << " Available result fields include: "
            << props.objectResultEntries(foName_)
            << exit(FatalError);
    }

    Type value = props.template getObjectResult<Type>(foName_, foResultName_);

    DebugInfo
        << "    Using " << foName_ << " function object value: " << value
        << endl;

    return value;
}


template<class Type>
inline Type Foam::Function1Types::FunctionObjectValue<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return pTraits<Type>::zero;
}


// ************************************************************************* //

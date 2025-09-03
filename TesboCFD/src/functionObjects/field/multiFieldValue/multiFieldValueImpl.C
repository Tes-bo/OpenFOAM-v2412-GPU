/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "FlatOutput.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Return magnitude of a given value unless the given Type is scalar
template<class Type>
static inline scalar MagOp(const Type& val)
{
    return Foam::mag(val);
}

//- Return self if a given Type is scalar
template<>
inline scalar MagOp<scalar>(const scalar& val)
{
    return val;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldValues::multiFieldValue::applyOperation
(
    const word& resultType,
    const wordList& foNames,
    const wordList& entryNames
)
{
    if (pTraits<Type>::typeName != resultType)
    {
        return false;
    }

    Type result = Zero;

    Field<Type> values(foNames.size());
    forAll(values, i)
    {
        values[i] = this->getObjectResult<Type>(foNames[i], entryNames[i]);
    }

    const word& opName = operationTypeNames_[operation_];

    switch (operation_)
    {
        case opSum:
        case opAdd:
        {
            result = sum(values);
            break;
        }
        case opSubtract:
        {
            result = values[0];
            for (label i = 1; i < values.size(); ++i)
            {
                result -= values[i];
            }
            break;
        }
        case opDivide:
        {
            result = values[0];
            for (label i = 1; i < values.size(); ++i)
            {
                result /= stabilise(MagOp(values[i]), SMALL);
            }
            break;
        }
        case opCmptDivide:
        {
            result = values[0];
            for (label i = 1; i < values.size(); ++i)
            {
                result = cmptDivide(result, stabilise(values[i], SMALL));
            }
            break;
        }
        case opMin:
        {
            result = min(values);
            break;
        }
        case opMax:
        {
            result = max(values);
            break;
        }
        case opAverage:
        {
            result = average(values);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unable to process operation "
                << operationTypeNames_[operation_]
                << abort(FatalError);
        }
    }

    OStringStream os;
    os << opName << flatOutput(entryNames, FlatOutput::ParenComma{});
    const word resultName(os.str());
    Log << "    " << resultName << " = " << result << endl;

    this->file()<< tab << result;

    // Write state/results information
    this->setResult(resultName, result);

    return true;
}


// ************************************************************************* //

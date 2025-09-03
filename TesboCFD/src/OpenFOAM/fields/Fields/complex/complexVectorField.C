/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "complexVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<complexVector>, complexVectorList);
    addCompoundToRunTimeSelectionTable(List<complexVector>, complexVectorList);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::zip
(
    complexVectorField& result,
    const UList<vector>& realValues,
    const UList<vector>& imagValues
)
{
    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != realValues.size() || len != imagValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << realValues.size() << ' ' << imagValues.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = Foam::zip(realValues[i], imagValues[i]);
    }
}


void Foam::zip
(
    complexVectorField& result,
    const UList<vector>& realValues,
    const vector& imagValue
)
{
    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != realValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " != "
            << realValues.size() << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = Foam::zip(realValues[i], imagValue);
    }
}


void Foam::zip
(
    complexVectorField& result,
    const vector& realValue,
    const UList<vector>& imagValues
)
{
    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != imagValues.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " != "
            << imagValues.size() << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = Foam::zip(realValue, imagValues[i]);
    }
}


Foam::complexVectorField Foam::ComplexField
(
    const UList<vector>& realValues,
    const UList<vector>& imagValues
)
{
    complexVectorField result(realValues.size());

    Foam::zip(result, realValues, imagValues);

    return result;
}


Foam::complexVectorField Foam::ComplexField
(
    const UList<vector>& realValues,
    const vector& imagValue
)
{
    complexVectorField result(realValues.size());

    Foam::zip(result, realValues, imagValue);

    return result;
}


Foam::complexVectorField Foam::ComplexField
(
    const vector& realValue,
    const UList<vector>& imagValues
)
{
    complexVectorField result(imagValues.size());

    Foam::zip(result, realValue, imagValues);

    return result;
}


Foam::vectorField Foam::ReImSum(const UList<complexVector>& cmplx)
{
    vectorField result(cmplx.size());

    std::transform
    (
        cmplx.cbegin(),
        cmplx.cend(),
        result.begin(),
        [](const complexVector& c)
        {
            return vector(c.x().cmptSum(), c.y().cmptSum(), c.z().cmptSum());
        }
    );

    return result;
}


Foam::vectorField Foam::Re(const UList<complexVector>& cmplx)
{
    vectorField result(cmplx.size());

    std::transform
    (
        cmplx.cbegin(),
        cmplx.cend(),
        result.begin(),
        [](const complexVector& c)
        {
            return vector(c.x().real(), c.y().real(), c.z().real());
        }
    );

    return result;
}


Foam::vectorField Foam::Im(const UList<complexVector>& cmplx)
{
    vectorField result(cmplx.size());

    std::transform
    (
        cmplx.cbegin(),
        cmplx.cend(),
        result.begin(),
        [](const complexVector& c)
        {
            return vector(c.x().imag(), c.y().imag(), c.z().imag());
        }
    );

    return result;
}


Foam::complexVectorField Foam::operator^
(
    const UList<vector>& vec,
    const UList<complexVector>& cmplx
)
{
    const label len = cmplx.size();

    #ifdef FULLDEBUG
    if (len != vec.size())
    {
        FatalErrorInFunction
            << "Parameter sizes do not match: " << vec.size()
            << " != " << len << nl
            << abort(FatalError);
    }
    #endif

    complexVectorField result(len);

    for (label i=0; i < len; ++i)
    {
        result[i] = (vec[i] ^ cmplx[i]);
    }

    return result;
}


// ************************************************************************* //

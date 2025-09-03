/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "vectorField.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt>
void Foam::zip
(
    Field<Vector<Cmpt>>& result,
    const UList<Cmpt>& x,
    const UList<Cmpt>& y,
    const UList<Cmpt>& z
)
{
    typedef Vector<Cmpt> value_type;

    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != x.size() || len != y.size() || len != z.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << x.size() << ' ' << y.size() << ' ' << z.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = value_type(x[i], y[i], z[i]);
    }
}


template<class Cmpt>
Foam::tmp<Foam::Field<Foam::Vector<Cmpt>>> Foam::zip
(
    const Field<Cmpt>& x,
    const Field<Cmpt>& y,
    const Field<Cmpt>& z
)
{
    auto tresult = tmp<Field<Vector<Cmpt>>>::New(x.size());

    Foam::zip(tresult.ref(), x, y, z);

    return tresult;
}


template<class Cmpt>
void Foam::unzip
(
    const UList<Vector<Cmpt>>& input,
    Field<Cmpt>& x,
    Field<Cmpt>& y,
    Field<Cmpt>& z
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if (len != x.size() || len != y.size() || len != z.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << x.size() << ' ' << y.size() << ' ' << z.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        x[i] = input[i].x();
        y[i] = input[i].y();
        z[i] = input[i].z();
    }
}


// ************************************************************************* //

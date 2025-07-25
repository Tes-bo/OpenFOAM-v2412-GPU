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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt>
void Foam::zip
(
    Field<SphericalTensor<Cmpt>>& result,
    const UList<Cmpt>& ii
)
{
    typedef SphericalTensor<Cmpt> value_type;

    const label len = result.size();

    #ifdef FULLDEBUG
    if (len != ii.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << ii.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = value_type(ii[i]);
    }
}


template<class Cmpt>
void Foam::unzip
(
    const UList<SphericalTensor<Cmpt>>& input,
    Field<Cmpt>& ii
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if (len != ii.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << ii.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        ii[i] = input[i].ii();
    }
}


template<class Cmpt>
Foam::tmp<Foam::Field<Foam::SphericalTensor<Cmpt>>>
Foam::zip
(
    const Field<Cmpt>& ii
)
{
    auto tresult = tmp<Field<SphericalTensor<Cmpt>>>::New(ii.size());

    Foam::zip(tresult.ref(), ii);

    return tresult;
}


// ************************************************************************* //

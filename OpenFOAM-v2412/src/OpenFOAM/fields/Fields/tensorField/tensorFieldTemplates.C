/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Cmpt>
void Foam::zip
(
    Field<Tensor<Cmpt>>& result,
    const UList<Cmpt>& xx, const UList<Cmpt>& xy, const UList<Cmpt>& xz,
    const UList<Cmpt>& yx, const UList<Cmpt>& yy, const UList<Cmpt>& yz,
    const UList<Cmpt>& zx, const UList<Cmpt>& zy, const UList<Cmpt>& zz
)
{
    typedef Tensor<Cmpt> value_type;

    const label len = result.size();

    #ifdef FULLDEBUG
    if
    (
        len != xx.size() || len != xy.size() || len != xz.size()
     || len != yx.size() || len != yy.size() || len != yz.size()
     || len != zx.size() || len != zy.size() || len != zz.size()
    )
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << xx.size() << ' ' << xy.size() << ' ' << xz.size() << ' '
            << yx.size() << ' ' << yy.size() << ' ' << yz.size() << ' '
            << zx.size() << ' ' << zy.size() << ' ' << zz.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = value_type
        (
            xx[i], xy[i], xz[i],
            yx[i], yy[i], yz[i],
            zx[i], zy[i], zz[i]
        );
    }
}


template<class Cmpt>
void Foam::unzip
(
    const UList<Tensor<Cmpt>>& input,
    Field<Cmpt>& xx, Field<Cmpt>& xy, Field<Cmpt>& xz,
    Field<Cmpt>& yx, Field<Cmpt>& yy, Field<Cmpt>& yz,
    Field<Cmpt>& zx, Field<Cmpt>& zy, Field<Cmpt>& zz
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if
    (
        len != xx.size() || len != xy.size() || len != xz.size()
     || len != yx.size() || len != yy.size() || len != yz.size()
     || len != zx.size() || len != zy.size() || len != zz.size()
    )
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << xx.size() << ' ' << xy.size() << ' ' << xz.size() << ' '
            << yx.size() << ' ' << yy.size() << ' ' << yz.size() << ' '
            << zx.size() << ' ' << zy.size() << ' ' << zz.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        xx[i] = input[i].xx(); xy[i] = input[i].xy(); xz[i] = input[i].xz();
        yx[i] = input[i].yx(); yy[i] = input[i].yy(); yz[i] = input[i].yz();
        zx[i] = input[i].zx(); zy[i] = input[i].zy(); zz[i] = input[i].zz();
    }
}


template<class Cmpt>
Foam::tmp<Foam::Field<Foam::Tensor<Cmpt>>>
Foam::zip
(
    const Field<Cmpt>& xx, const Field<Cmpt>& xy, const Field<Cmpt>& xz,
    const Field<Cmpt>& yx, const Field<Cmpt>& yy, const Field<Cmpt>& yz,
    const Field<Cmpt>& zx, const Field<Cmpt>& zy, const Field<Cmpt>& zz
)
{
    auto tresult = tmp<Field<Tensor<Cmpt>>>::New(xx.size());

    Foam::zip(tresult.ref(), xx, xy, xz, yx, yy, yz, zx, zy, zz);

    return tresult;
}


template<class Cmpt>
void Foam::zipRows
(
    Field<Tensor<Cmpt>>& result,
    const UList<Vector<Cmpt>>& x,
    const UList<Vector<Cmpt>>& y,
    const UList<Vector<Cmpt>>& z
)
{
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
        result[i].rows(x[i], y[i], z[i]);
    }
}


template<class Cmpt>
void Foam::zipCols
(
    Field<Tensor<Cmpt>>& result,
    const UList<Vector<Cmpt>>& x,
    const UList<Vector<Cmpt>>& y,
    const UList<Vector<Cmpt>>& z
)
{
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
        result[i].cols(x[i], y[i], z[i]);
    }
}


template<class Cmpt>
void Foam::unzipRows
(
    const UList<Tensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& x,
    Field<Vector<Cmpt>>& y,
    Field<Vector<Cmpt>>& z
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


template<class Cmpt>
void Foam::unzipCols
(
    const UList<Tensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& x,
    Field<Vector<Cmpt>>& y,
    Field<Vector<Cmpt>>& z
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
        x[i] = input[i].cx();
        y[i] = input[i].cy();
        z[i] = input[i].cz();
    }
}


template<class Cmpt>
void Foam::unzipRow
(
    const UList<Tensor<Cmpt>>& input,
    const direction idx,
    Field<Vector<Cmpt>>& result
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if (len != result.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << result.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    switch (idx)
    {
        case vector::components::X :
        {
            for (label i=0; i < len; ++i)
            {
                result[i] = input[i].x();
            }
        }
        break;

        case vector::components::Y :
        {
            for (label i=0; i < len; ++i)
            {
                result[i] = input[i].y();
            }
        }
        break;

        case vector::components::Z :
        {
            for (label i=0; i < len; ++i)
            {
                result[i] = input[i].z();
            }
        }
        break;
    }
}


template<class Cmpt>
void Foam::unzipCol
(
    const UList<Tensor<Cmpt>>& input,
    const direction idx,
    Field<Vector<Cmpt>>& result
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if (len != result.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << result.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    switch (idx)
    {
        case vector::components::X :
        {
            for (label i=0; i < len; ++i)
            {
                result[i] = input[i].cx();
            }
        }
        break;

        case vector::components::Y :
        {
            for (label i=0; i < len; ++i)
            {
                result[i] = input[i].cy();
            }
        }
        break;

        case vector::components::Z :
        {
            for (label i=0; i < len; ++i)
            {
                result[i] = input[i].cz();
            }
        }
        break;
    }
}


template<class Cmpt>
void Foam::unzipDiag
(
    const UList<Tensor<Cmpt>>& input,
    Field<Vector<Cmpt>>& result
)
{
    const label len = input.size();

    #ifdef FULLDEBUG
    if (len != result.size())
    {
        FatalErrorInFunction
            << "Components sizes do not match: " << len << " ("
            << result.size() << ')'
            << nl
            << abort(FatalError);
    }
    #endif

    for (label i=0; i < len; ++i)
    {
        result[i] = input[i].diag();
    }
}


template<class Cmpt>
Foam::tmp<Foam::Field<Foam::Vector<Cmpt>>>
Foam::unzipRow
(
    const Field<Tensor<Cmpt>>& input,
    const direction idx
)
{
    auto tresult = tmp<Field<Vector<Cmpt>>>::New(input.size());

    Foam::unzipRow(input, idx, tresult.ref());

    return tresult;
}


template<class Cmpt>
Foam::tmp<Foam::Field<Foam::Vector<Cmpt>>>
Foam::unzipCol
(
    const Field<Tensor<Cmpt>>& input,
    const direction idx
)
{
    auto tresult = tmp<Field<Vector<Cmpt>>>::New(input.size());

    Foam::unzipCol(input, idx, tresult.ref());

    return tresult;
}


template<class Cmpt>
Foam::tmp<Foam::Field<Foam::Vector<Cmpt>>>
Foam::unzipDiag
(
    const Field<Tensor<Cmpt>>& input
)
{
    auto tresult = tmp<Field<Vector<Cmpt>>>::New(input.size());

    Foam::unzipDiag(input, tresult.ref());

    return tresult;
}


// ************************************************************************* //

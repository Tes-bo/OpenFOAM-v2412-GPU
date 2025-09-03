/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "DimensionedScalarField.H"

#define TEMPLATE template<class GeoMesh>
#include "DimensionedFieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> stabilise
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const dimensioned<scalar>& ds
)
{
    auto tres =
        DimensionedField<scalar, GeoMesh>::New
        (
            "stabilise(" + dsf.name() + ',' + ds.name() + ')',
            dsf.mesh(),
            dsf.dimensions() + ds.dimensions()
        );

    stabilise(tres.ref().field(), dsf.field(), ds.value());

    return tres;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> stabilise
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf,
    const dimensioned<scalar>& ds
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh>> tres = New
    (
        tdsf,
        "stabilise(" + dsf.name() + ',' + ds.name() + ')',
        dsf.dimensions() + ds.dimensions()
    );

    stabilise(tres.ref().field(), dsf.field(), ds.value());

    tdsf.clear();

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, '+', add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, '-', subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, '*', multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, '|', divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, '|', divide)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const DimensionedField<scalar, GeoMesh>& f1,
    const DimensionedField<scalar, GeoMesh>& f2
)
{
    if
    (
        dimensionSet::checking()
     && (!f1.dimensions().dimensionless() || !f2.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl;

        if (!f1.dimensions().dimensionless())
        {
            FatalError
                << "    Base field dimensions: " << f1.dimensions() << nl;
        }
        if (!f2.dimensions().dimensionless())
        {
            FatalError
                << "    Exponent field dimensions: " << f2.dimensions() << nl;
        }
        FatalError << exit(FatalError);
    }

    auto tresult =
        DimensionedField<scalar, GeoMesh>::New
        (
            "pow(" + f1.name() + ',' + f2.name() + ')',
            f1.mesh(),
            dimless
        );

    pow(tresult.ref().field(), f1.field(), f2.field());

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tf1,
    const DimensionedField<scalar, GeoMesh>& f2
)
{
    const auto& f1 = tf1();

    if
    (
        dimensionSet::checking()
     && (!f1.dimensions().dimensionless() || !f2.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl;

        if (!f1.dimensions().dimensionless())
        {
            FatalError
                << "    Base field dimensions: " << f1.dimensions() << nl;
        }
        if (!f2.dimensions().dimensionless())
        {
            FatalError
                << "    Exponent field dimensions: " << f2.dimensions() << nl;
        }
        FatalError << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh>> tresult = New
    (
        tf1,
        "pow(" + f1.name() + ',' + f2.name() + ')',
        dimless
    );

    pow(tresult.ref().field(), f1.field(), f2.field());

    tf1.clear();

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const DimensionedField<scalar, GeoMesh>& f1,
    const tmp<DimensionedField<scalar, GeoMesh>>& tf2
)
{
    const auto& f2 = tf2();

    if
    (
        dimensionSet::checking()
     && (!f1.dimensions().dimensionless() || !f2.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl;

        if (!f1.dimensions().dimensionless())
        {
            FatalError
                << "    Base field dimensions: " << f1.dimensions() << nl;
        }
        if (!f2.dimensions().dimensionless())
        {
            FatalError
                << "    Exponent field dimensions: " << f2.dimensions() << nl;
        }
        FatalError << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh>> tresult = New
    (
        tf2,
        "pow(" + f1.name() + ',' + f2.name() + ')',
        dimless
    );

    pow(tresult.ref().field(), f1.field(), f2.field());

    tf2.clear();

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tf1,
    const tmp<DimensionedField<scalar, GeoMesh>>& tf2
)
{
    const auto& f1 = tf1();
    const auto& f2 = tf2();

    if
    (
        dimensionSet::checking()
     && (!f1.dimensions().dimensionless() || !f2.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl;

        if (!f1.dimensions().dimensionless())
        {
            FatalError
                << "    Base field dimensions: " << f1.dimensions() << nl;
        }
        if (!f2.dimensions().dimensionless())
        {
            FatalError
                << "    Exponent field dimensions: " << f2.dimensions() << nl;
        }
        FatalError << exit(FatalError);
    }

    auto tresult =
        reuseTmpTmpDimensionedField<scalar, scalar, scalar, scalar, GeoMesh>::
        New
        (
            tf1,
            tf2,
            "pow(" + f1.name() + ',' + f2.name() + ')',
            dimless
        );

    pow(tresult.ref().field(), f1.field(), f2.field());

    tf1.clear();
    tf2.clear();

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const DimensionedField<scalar, GeoMesh>& f1,
    const dimensionedScalar& ds
)
{
    if
    (
        dimensionSet::checking()
     && (!ds.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl
            << "    Exponent dimensions: " << ds.dimensions() << nl
            << exit(FatalError);
    }

    auto tresult =
        DimensionedField<scalar, GeoMesh>::New
        (
            "pow(" + f1.name() + ',' + ds.name() + ')',
            f1.mesh(),
            pow(f1.dimensions(), ds)
        );

    pow(tresult.ref().field(), f1.field(), ds.value());

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tf1,
    const dimensionedScalar& ds
)
{
    const auto& f1 = tf1();

    if
    (
        dimensionSet::checking()
     && (!ds.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl
            << "    Exponent dimensions: " << ds.dimensions() << nl
            << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh>> tresult = New
    (
        tf1,
        "pow(" + f1.name() + ',' + ds.name() + ')',
        pow(f1.dimensions(), ds)
    );

    pow(tresult.ref().field(), f1.field(), ds.value());

    tf1.clear();

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const scalar& s
)
{
    return pow(dsf, dimensionedScalar(s));
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf,
    const scalar& s
)
{
    return pow(tdsf, dimensionedScalar(s));
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const dimensionedScalar& ds,
    const DimensionedField<scalar, GeoMesh>& f2
)
{
    if
    (
        dimensionSet::checking()
     && (!ds.dimensions().dimensionless() || !f2.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl;

        if (!ds.dimensions().dimensionless())
        {
            FatalError
                << "    Base scalar dimensions: " << ds.dimensions() << nl;
        }
        if (!f2.dimensions().dimensionless())
        {
            FatalError
                << "    Exponent field dimensions: " << f2.dimensions() << nl;
        }
        FatalError << exit(FatalError);
    }

    auto tresult =
        DimensionedField<scalar, GeoMesh>::New
        (
            "pow(" + ds.name() + ',' + f2.name() + ')',
            f2.mesh(),
            dimless
        );

    pow(tresult.ref().field(), ds.value(), f2.field());

    return tresult;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const dimensionedScalar& ds,
    const tmp<DimensionedField<scalar, GeoMesh>>& tf2
)
{
    const auto& f2 = tf2();

    if
    (
        dimensionSet::checking()
     && (!ds.dimensions().dimensionless() || !f2.dimensions().dimensionless())
    )
    {
        FatalErrorInFunction
            << "pow() expects dimensionless parameters, but found" << nl;

        if (!ds.dimensions().dimensionless())
        {
            FatalError
                << "    Base scalar dimensions: " << ds.dimensions() << nl;
        }
        if (!f2.dimensions().dimensionless())
        {
            FatalError
                << "    Exponent field dimensions: " << f2.dimensions() << nl;
        }

        FatalError << exit(FatalError);
    }

    tmp<DimensionedField<scalar, GeoMesh>> tresult = New
    (
        tf2,
        "pow(" + ds.name() + ',' + f2.name() + ')',
        dimless
    );

    pow(tresult.ref().field(), ds.value(), f2.field());

    tf2.clear();

    return tresult;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const scalar& s,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    return pow(dimensionedScalar(s), dsf);
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> pow
(
    const scalar& s,
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf
)
{
    return pow(dimensionedScalar(s), tdsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf1,
    const DimensionedField<scalar, GeoMesh>& dsf2
)
{
    auto tres =
        DimensionedField<scalar, GeoMesh>::New
        (
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            dsf1.mesh(),
            atan2(dsf1.dimensions(), dsf2.dimensions())
        );

    atan2(tres.ref().field(), dsf1.field(), dsf2.field());

    return tres;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf1,
    const DimensionedField<scalar, GeoMesh>& dsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf1 = tdsf1();

    tmp<DimensionedField<scalar, GeoMesh>> tres = New
    (
        tdsf1,
        "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
        atan2(dsf1.dimensions(), dsf2.dimensions())
    );

    atan2(tres.ref().field(), dsf1.field(), dsf2.field());

    tdsf1.clear();

    return tres;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf1,
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf2 = tdsf2();

    tmp<DimensionedField<scalar, GeoMesh>> tres = New
    (
        tdsf2,
        "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
        atan2(dsf1.dimensions(), dsf2.dimensions())
    );

    atan2(tres.ref().field(), dsf1.field(), dsf2.field());

    tdsf2.clear();

    return tres;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf1,
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf2
)
{
    const DimensionedField<scalar, GeoMesh>& dsf1 = tdsf1();
    const DimensionedField<scalar, GeoMesh>& dsf2 = tdsf2();

    auto tres =
        reuseTmpTmpDimensionedField<scalar, scalar, scalar, scalar, GeoMesh>::
        New
        (
            tdsf1,
            tdsf2,
            "atan2(" + dsf1.name() + ',' + dsf2.name() + ')',
            atan2(dsf1.dimensions(), dsf2.dimensions())
        );

    atan2(tres.ref().field(), dsf1.field(), dsf2.field());

    tdsf1.clear();
    tdsf2.clear();

    return tres;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const dimensionedScalar& ds
)
{
    auto tres =
        DimensionedField<scalar, GeoMesh>::New
        (
            "atan2(" + dsf.name() + ',' + ds.name() + ')',
            dsf.mesh(),
            atan2(dsf.dimensions(), ds)
        );

    atan2(tres.ref().field(), dsf.field(), ds.value());

    return tres;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf,
    const dimensionedScalar& ds
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh>> tres = New
    (
        tdsf,
        "atan2(" + dsf.name() + ',' + ds.name() + ')',
        atan2(dsf.dimensions(), ds)
    );

    atan2(tres.ref().field(), dsf.field(), ds.value());

    tdsf.clear();

    return tres;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const DimensionedField<scalar, GeoMesh>& dsf,
    const scalar& s
)
{
    return atan2(dsf, dimensionedScalar(s));
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf,
    const scalar& s
)
{
    return atan2(tdsf, dimensionedScalar(s));
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const dimensionedScalar& ds,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    auto tres =
        DimensionedField<scalar, GeoMesh>::New
        (
            "atan2(" + ds.name() + ',' + dsf.name() + ')',
            dsf.mesh(),
            atan2(ds, dsf.dimensions())
        );

    atan2(tres.ref().field(), ds.value(), dsf.field());

    return tres;
}


template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const dimensionedScalar& ds,
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf
)
{
    const DimensionedField<scalar, GeoMesh>& dsf = tdsf();

    tmp<DimensionedField<scalar, GeoMesh>> tres = New
    (
        tdsf,
        "atan2(" + ds.name() + ',' + dsf.name() + ')',
        atan2(ds, dsf.dimensions())
    );

    atan2(tres.ref().field(), ds.value(), dsf.field());

    tdsf.clear();

    return tres;
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const scalar& s,
    const DimensionedField<scalar, GeoMesh>& dsf
)
{
    return atan2(dimensionedScalar(s), dsf);
}

template<class GeoMesh>
tmp<DimensionedField<scalar, GeoMesh>> atan2
(
    const scalar& s,
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf
)
{
    return atan2(dimensionedScalar(s), tdsf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, scalar, pow3, pow3)
UNARY_FUNCTION(scalar, scalar, pow4, pow4)
UNARY_FUNCTION(scalar, scalar, pow5, pow5)
UNARY_FUNCTION(scalar, scalar, pow6, pow6)
UNARY_FUNCTION(scalar, scalar, pow025, pow025)
UNARY_FUNCTION(scalar, scalar, sqrt, sqrt)
UNARY_FUNCTION(scalar, scalar, cbrt, cbrt)
UNARY_FUNCTION(scalar, scalar, sign, sign)
UNARY_FUNCTION(scalar, scalar, pos, pos)
UNARY_FUNCTION(scalar, scalar, pos0, pos0)
UNARY_FUNCTION(scalar, scalar, neg, neg)
UNARY_FUNCTION(scalar, scalar, neg0, neg0)
UNARY_FUNCTION(scalar, scalar, posPart, posPart)
UNARY_FUNCTION(scalar, scalar, negPart, negPart)

UNARY_FUNCTION(scalar, scalar, exp, trans)
UNARY_FUNCTION(scalar, scalar, log, trans)
UNARY_FUNCTION(scalar, scalar, log10, trans)
UNARY_FUNCTION(scalar, scalar, sin, trans)
UNARY_FUNCTION(scalar, scalar, cos, trans)
UNARY_FUNCTION(scalar, scalar, tan, trans)
UNARY_FUNCTION(scalar, scalar, asin, trans)
UNARY_FUNCTION(scalar, scalar, acos, trans)
UNARY_FUNCTION(scalar, scalar, atan, trans)
UNARY_FUNCTION(scalar, scalar, sinh, trans)
UNARY_FUNCTION(scalar, scalar, cosh, trans)
UNARY_FUNCTION(scalar, scalar, tanh, trans)
UNARY_FUNCTION(scalar, scalar, asinh, trans)
UNARY_FUNCTION(scalar, scalar, acosh, trans)
UNARY_FUNCTION(scalar, scalar, atanh, trans)
UNARY_FUNCTION(scalar, scalar, erf, trans)
UNARY_FUNCTION(scalar, scalar, erfc, trans)
UNARY_FUNCTION(scalar, scalar, lgamma, trans)
UNARY_FUNCTION(scalar, scalar, j0, trans)
UNARY_FUNCTION(scalar, scalar, j1, trans)
UNARY_FUNCTION(scalar, scalar, y0, trans)
UNARY_FUNCTION(scalar, scalar, y1, trans)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define BesselFunc(func)                                                       \
                                                                               \
template<class GeoMesh>                                                        \
tmp<DimensionedField<scalar, GeoMesh>> func                                    \
(                                                                              \
    const int n,                                                               \
    const DimensionedField<scalar, GeoMesh>& dsf                               \
)                                                                              \
{                                                                              \
    if (dimensionSet::checking() && !dsf.dimensions().dimensionless())         \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "Field is not dimensionless: " << dsf.dimensions() << nl        \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    auto tres =                                                                \
        DimensionedField<scalar, GeoMesh>::New                                 \
        (                                                                      \
            #func "(" + name(n) + ',' + dsf.name() + ')',                      \
            dsf.mesh(),                                                        \
            dimless                                                            \
        );                                                                     \
                                                                               \
    func(tres.ref().field(), n, dsf.field());                                  \
                                                                               \
    return tres;                                                               \
}                                                                              \
                                                                               \
                                                                               \
template<class GeoMesh>                                                        \
tmp<DimensionedField<scalar, GeoMesh>> func                                    \
(                                                                              \
    const int n,                                                               \
    const tmp<DimensionedField<scalar, GeoMesh>>& tdsf                         \
)                                                                              \
{                                                                              \
    const auto& dsf = tdsf();                                                  \
                                                                               \
    if (dimensionSet::checking() && !dsf.dimensions().dimensionless())         \
    {                                                                          \
        FatalErrorInFunction                                                   \
            << "Field is not dimensionless: " << dsf.dimensions() << nl        \
            << abort(FatalError);                                              \
    }                                                                          \
                                                                               \
    tmp<DimensionedField<scalar, GeoMesh>> tres                                \
    (                                                                          \
        New                                                                    \
        (                                                                      \
            tdsf,                                                              \
            #func "(" + name(n) + ',' + dsf.name() + ')',                      \
            dimless                                                            \
        )                                                                      \
    );                                                                         \
                                                                               \
    func(tres.ref().field(), n, dsf.field());                                  \
                                                                               \
    tdsf.clear();                                                              \
                                                                               \
    return tres;                                                               \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "facDiv.H"
#include "faMesh.H"
#include "facEdgeIntegrate.H"
#include "faDivScheme.H"
#include "faConvectionScheme.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const GeometricField<Type, faePatchField, edgeMesh>& ssf
)
{
    const areaVectorField& n = ssf.mesh().faceAreaNormals();

    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv =
        fac::edgeIntegrate(ssf);

    GeometricField<Type, faPatchField, areaMesh>& Div = tDiv.ref();

    Div.primitiveFieldRef() =
        transform(tensor::I - sqr(n), Div.internalField());
    Div.correctBoundaryConditions();

    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<GeometricField<Type, faePatchField, edgeMesh>>& tssf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv(fac::div(tssf()));
    tssf.clear();
    return tDiv;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    const areaVectorField& n = vf.mesh().faceAreaNormals();

    tmp
    <
        GeometricField
        <
            typename innerProduct<vector, Type>::type,
            faPatchField,
            areaMesh
        >
    > tDiv
    (
        fa::divScheme<Type>::New
        (
            vf.mesh(), vf.mesh().divScheme(name)
        ).ref().facDiv(vf)
    );
    GeometricField
    <
        typename innerProduct<vector, Type>::type,
        faPatchField,
        areaMesh
    >& Div = tDiv.ref();

    Div.primitiveFieldRef() =
        transform(tensor::I - sqr(n), Div.internalField());

    Div.correctBoundaryConditions();

    return tDiv;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvvf,
    const word& name
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh>> tDiv
    (
        fac::div(tvvf(), name)
    );
    tvvf.clear();
    return tDiv;
}

template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::div(vf, "div("+vf.name()+')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename innerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
div
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvvf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    tmp<GeometricField<DivType, faPatchField, areaMesh>> tDiv
    (
        fac::div(tvvf())
    );

    tvvf.clear();
    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    const areaVectorField& n = vf.mesh().faceAreaNormals();

    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fa::convectionScheme<Type>::New
        (
            vf.mesh(),
            flux,
            vf.mesh().divScheme(name)
        ).ref().facDiv(flux, vf)
    );
    GeometricField<Type, faPatchField, areaMesh>& Div = tDiv.ref();

    Div.primitiveFieldRef() = transform(tensor::I - sqr(n), Div.internalField());
    Div.correctBoundaryConditions();

    return tDiv;

}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fac::div(tflux(), vf, name)
    );
    tflux.clear();
    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fac::div(flux, tvf(), name)
    );
    tvf.clear();
    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fac::div(tflux(), tvf(), name)
    );
    tflux.clear();
    tvf.clear();
    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::div
    (
        flux, vf, "div("+flux.name()+','+vf.name()+')'
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fac::div(tflux(), vf)
    );
    tflux.clear();
    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const edgeScalarField& flux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fac::div(flux, tvf())
    );
    tvf.clear();
    return tDiv;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
div
(
    const tmp<edgeScalarField>& tflux,
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh>> tDiv
    (
        fac::div(tflux(), tvf())
    );
    tflux.clear();
    tvf.clear();
    return tDiv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

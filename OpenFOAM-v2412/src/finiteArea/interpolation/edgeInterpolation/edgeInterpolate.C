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

#include "edgeInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::edgeInterpolationScheme<Type>> Foam::fac::scheme
(
    const edgeScalarField& faceFlux,
    Istream& streamData
)
{
    return edgeInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        streamData
    );
}


template<class Type>
Foam::tmp<Foam::edgeInterpolationScheme<Type>> Foam::fac::scheme
(
    const edgeScalarField& faceFlux,
    const word& name
)
{
    return edgeInterpolationScheme<Type>::New
    (
        faceFlux.mesh(),
        faceFlux,
        faceFlux.mesh().interpolationScheme(name)
    );
}


template<class Type>
Foam::tmp<Foam::edgeInterpolationScheme<Type>> Foam::fac::scheme
(
    const faMesh& mesh,
    Istream& streamData
)
{
    return edgeInterpolationScheme<Type>::New
    (
        mesh,
        streamData
    );
}


template<class Type>
Foam::tmp<Foam::edgeInterpolationScheme<Type>> Foam::fac::scheme
(
    const faMesh& mesh,
    const word& name
)
{
    return edgeInterpolationScheme<Type>::New
    (
        mesh,
        mesh.interpolationScheme(name)
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const edgeScalarField& faceFlux,
    Istream& schemeData
)
{
    #ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << endl;
    }
    #endif

    return scheme<Type>(faceFlux, schemeData)().interpolate(vf);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const edgeScalarField& faceFlux,
    const word& name
)
{
    #ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using " << name
            << endl;
    }
    #endif

    return scheme<Type>(faceFlux, name)().interpolate(vf);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const edgeScalarField& faceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        interpolate(tvf(), faceFlux, name);

    tvf.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const tmp<edgeScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        interpolate(vf, tFaceFlux(), name);

    tFaceFlux.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const tmp<edgeScalarField>& tFaceFlux,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        interpolate(tvf(), tFaceFlux(), name);

    tvf.clear();
    tFaceFlux.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    Istream& schemeData
)
{
    #ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << endl;
    }
    #endif

    return scheme<Type>(vf.mesh(), schemeData)().interpolate(vf);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    #ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using " << name
            << endl;
    }
    #endif

    return scheme<Type>(vf.mesh(), name)().interpolate(vf);
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        interpolate(tvf(), name);

    tvf.clear();

    return tsf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    #ifdef DEBUGInterpolations
    if (edgeInterpolation::debug)
    {
        InfoInFunction
            << "interpolating GeometricField<Type, faPatchField, areaMesh> "
            << "using run-time selected scheme"
            << endl;
    }
    #endif

    return interpolate(vf, "interpolate(" + vf.name() + ')');
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh>>
Foam::fac::interpolate
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& tvf
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh>> tsf =
        interpolate(tvf());
    tvf.clear();

    return tsf;
}


// ************************************************************************* //

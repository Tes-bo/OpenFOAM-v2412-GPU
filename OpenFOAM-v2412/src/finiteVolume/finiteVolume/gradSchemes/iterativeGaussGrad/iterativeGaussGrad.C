/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "iterativeGaussGrad.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp
<
    Foam::GeometricField
    <
        typename Foam::outerProduct<Foam::vector, Type>::type,
        Foam::fvPatchField,
        Foam::volMesh
    >
>
Foam::fv::iterativeGaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;
    typedef GeometricField<GradType, fvsPatchField, surfaceMesh>
        GradSurfFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfFieldType;

    tmp<SurfFieldType> tssf = linearInterpolate(vsf);
    const SurfFieldType& ssf = tssf.cref();

    tmp<GradFieldType> tgGrad = fv::gaussGrad<Type>::gradf(ssf, name);
    GradFieldType& gGrad = tgGrad.ref();

    const skewCorrectionVectors& skv = skewCorrectionVectors::New(vsf.mesh());

    scalar relax = 1;
    const bool useRelax =
        vsf.mesh().relaxField("grad(" + vsf.name() + ")", relax);

    for (label i = 0; i < nIter_; ++i)
    {
        tmp<GradSurfFieldType> tsgGrad = linearInterpolate(gGrad);

        tmp<SurfFieldType> tcorr = skv() & tsgGrad;

        tcorr.ref().dimensions().reset(vsf.dimensions());

        if (useRelax)
        {
            // relax*prediction + (1-relax)*old
            gGrad *= (1.0 - relax);
            gGrad += relax*fv::gaussGrad<Type>::gradf(tcorr + ssf, name);
        }
        else
        {
            gGrad = fv::gaussGrad<Type>::gradf(tcorr + ssf, name);
        }
    }

    fv::gaussGrad<Type>::correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "deferredCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::deferredCorrection<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfCorr
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            "deferredCorrection::correction(" + vf.name() + ')',
            tbaseScheme_().interpolate(vf)
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr.ref();

    // Interpolate using the upwind weights to avoid circular reference to
    // [this] explicit correction
    sfCorr -= upwind<Type>::interpolate(vf, upwind<Type>::weights());
/*
    auto& sfCorrBf = sfCorr.boundaryFieldRef();
    for (auto& pf : sfCorrBf)
    {
        if (!pf.coupled())
        {
            pf = pTraits<Type>::zero;
        }
    }
*/
    return tsfCorr;
}


namespace Foam
{
    //makeSurfaceInterpolationScheme(deferredCorrection)
    makeSurfaceInterpolationTypeScheme(deferredCorrection, scalar)
    makeSurfaceInterpolationTypeScheme(deferredCorrection, vector)
}

// ************************************************************************* //

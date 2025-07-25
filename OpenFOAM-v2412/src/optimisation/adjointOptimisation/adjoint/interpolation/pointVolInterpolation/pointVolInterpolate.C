/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) Wikki Ltd
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

#include "pointVolInterpolation.H"
#include "volFields.H"
#include "pointFields.H"
#include "primitiveMesh.H"
#include "emptyFvPatch.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::pointVolInterpolation::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    DebugInFunction
        << "interpolating field from points to cells"
        << endl;

    const FieldField<Field, scalar>& weights = volWeights();
    const labelListList& cellPoints = vf.mesh().cellPoints();

    // Multiply pointField by weighting factor matrix to create volField
    forAll(cellPoints, cellI)
    {
        vf[cellI] = pTraits<Type>::zero;

        const labelList& curCellPoints = cellPoints[cellI];

        forAll(curCellPoints, cellPointI)
        {
            vf[cellI] +=
                weights[cellI][cellPointI]*pf[curCellPoints[cellPointI]];
        }
    }


    // Interpolate patch values: over-ride the internal values for the points
    // on the patch with the interpolated point values from the faces
    const fvBoundaryMesh& bm = vMesh().boundary();

    const PtrList<primitivePatchInterpolation>& pi = patchInterpolators();
    forAll(bm, patchI)
    {
        // If the patch is empty, skip it
        if (bm[patchI].type() != emptyFvPatch::typeName)
        {
            vf.boundaryFieldRef()[patchI] =
                pi[patchI].pointToFaceInterpolate
                (
                    pf.boundaryField()[patchI].patchInternalField()
                );
        }
    }

    vf.correctBoundaryConditions();

    DebugInFunction
        << "finished interpolating field from points to cells"
        << endl;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::pointVolInterpolation::interpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    // Construct tmp<pointField>
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "pointVolInterpolate(" + pf.name() + ')',
                pf.instance(),
                pf.db()
            ),
            vMesh(),
            pf.dimensions()
        )
    );

    // Perform interpolation
    interpolate(pf, tvf.ref());

    return tvf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::pointVolInterpolation::interpolate
(
    const tmp<GeometricField<Type, pointPatchField, pointMesh>>& tpf
) const
{
    // Construct tmp<volField>
    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf =
        interpolate(tpf());
    tpf.clear();
    return tvf;
}


// ************************************************************************* //

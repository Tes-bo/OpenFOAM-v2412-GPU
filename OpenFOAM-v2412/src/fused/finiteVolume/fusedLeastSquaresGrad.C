/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
    Copyright (C) 2024 M. Janssens
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

#include "fusedLeastSquaresGrad.H"
#include "leastSquaresVectors.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "extrapolatedCalculatedFvPatchField.H"
#include "fvcSurfaceOps.H"

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
Foam::fv::fusedLeastSquaresGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = vf.mesh();

    DebugPout<< "fusedLeastSquaresGrad<Type>::calcGrad on " << vf.name()
        << " with name " << name << endl;

    tmp<GradFieldType> tlsGrad
    (
        new GradFieldType
        (
            IOobject
            (
                name,
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>(vf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GradFieldType& lsGrad = tlsGrad.ref();

    // Get reference to least square vectors
    const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);
    const surfaceVectorField& ownLs = lsv.pVectors();
    const surfaceVectorField& neiLs = lsv.nVectors();

    const auto differenceOp = [&]
    (
        const vector& area,
        const Type& ownVal,
        const Type& neiVal
    ) -> Type
    {
        return neiVal - ownVal;
    };

    fvc::surfaceOp
    (
        vf,
        ownLs,
        neiLs,
        differenceOp,
        lsGrad
    );

    gaussGrad<Type>::correctBoundaryConditions(vf, lsGrad);

    return tlsGrad;
}


// ************************************************************************* //

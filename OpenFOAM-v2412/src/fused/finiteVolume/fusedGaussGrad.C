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

#include "fusedGaussGrad.H"
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
Foam::fv::fusedGaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = ssf.mesh();

    DebugPout<< "fusedGaussGrad<Type>::gradf on " << ssf.name()
        << " with name " << name << endl;

    tmp<GradFieldType> tgGrad
    (
        new GradFieldType
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>(ssf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GradFieldType& gGrad = tgGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();

    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    forAll(owner, facei)
    {
        const GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];
        }
    }

    igGrad /= mesh.V();

    gGrad.correctBoundaryConditions();

    return tgGrad;
}


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
Foam::fv::fusedGaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const fvMesh& mesh = vf.mesh();

    DebugPout<< "fusedGaussGrad<Type>::calcGrad on " << vf.name()
        << " with name " << name << endl;

    tmp<GradFieldType> tgGrad
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
    GradFieldType& gGrad = tgGrad.ref();

    if (this->tinterpScheme_().corrected())
    {
        const auto tfaceCorr(this->tinterpScheme_().correction(vf));
        auto& faceCorr = tfaceCorr();

        DebugPout<< "fusedGaussGrad<Type>::calcGrad corrected interpScheme "
            << this->tinterpScheme_().type() << endl;

        const auto interpolate = [&]
        (
            const vector& area,
            const scalar lambda,

            const Type& ownVal,
            const Type& neiVal,

            const Type& correction

        ) -> GradType
        {
            return area*((lambda*(ownVal - neiVal) + neiVal) + correction);
        };

        fvc::surfaceSum
        (
            this->tinterpScheme_().weights(vf),
            vf,
            faceCorr,
            interpolate,
            gGrad,
            false
        );
    }
    else
    {
        DebugPout<< "fusedGaussGrad<Type>::calcGrad uncorrected interpScheme "
            << this->tinterpScheme_().type() << endl;

        const auto interpolate = [&]
        (
            const vector& area,
            const scalar lambda,
            const Type& ownVal,
            const Type& neiVal
        ) -> GradType
        {
            return area*(lambda*(ownVal - neiVal) + neiVal);
        };

        fvc::surfaceSum
        (
            tinterpScheme_().weights(vf),
            vf,
            interpolate,
            gGrad,
            false
        );
    }

    gGrad.primitiveFieldRef() /= mesh.V();

    gGrad.correctBoundaryConditions();

    correctBoundaryConditions(vf, gGrad);

    return tgGrad;
}


template<class Type>
template<class GradType>
void Foam::fv::fusedGaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    GeometricField<GradType, fvPatchField, volMesh>& gGrad
)
{
    DebugPout<< "fusedGaussGrad<Type>::correctBoundaryConditions on "
        << vf.name() << " with gGrad " << gGrad.name() << endl;

    const fvMesh& mesh = vf.mesh();
    auto& gGradbf = gGrad.boundaryFieldRef();

    forAll(vf.boundaryField(), patchi)
    {
        if (!vf.boundaryField()[patchi].coupled())
        {
            const auto& pSf = mesh.Sf().boundaryField()[patchi];
            const auto tsnGrad(vf.boundaryField()[patchi].snGrad());
            const auto& snGrad = tsnGrad();
            auto& pgrad = gGradbf[patchi];

            forAll(pgrad, facei)
            {
                const vector n(pSf[facei]/mag(pSf[facei]));
                const Type uncorrectSnGrad(n & pgrad[facei]);
                pgrad[facei] += n*(snGrad[facei] - uncorrectSnGrad);
            }
        }
    }
}


// ************************************************************************* //

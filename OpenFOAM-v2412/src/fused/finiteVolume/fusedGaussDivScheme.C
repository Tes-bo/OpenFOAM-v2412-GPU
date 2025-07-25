/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvcSurfaceOps.H"
#include "fusedGaussDivScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
>
fusedGaussDivScheme<Type>::fvcDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef typename innerProduct<vector, Type>::type DivType;
    typedef GeometricField<DivType, fvPatchField, volMesh> DivFieldType;

    const fvMesh& mesh = vf.mesh();

    DebugPout<< "fusedGaussDivScheme<Type>::fvcDiv on " << vf.name()
         << endl;

    tmp<DivFieldType> tDiv
    (
        new DivFieldType
        (
            IOobject
            (
                "div(" + vf.name() + ')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<DivType>(vf.dimensions()/dimLength, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    DivFieldType& div = tDiv.ref();

    //fvcDiv(vf, this->tinterpScheme_().weights(vf), div);

    if (this->tinterpScheme_().corrected())
    {
        const auto tfaceCorr(this->tinterpScheme_().correction(vf));
        auto& faceCorr = tfaceCorr();

        const auto dotInterpolate = [&]
        (
            const vector& area,
            const scalar lambda,

            const Type& ownVal,
            const Type& neiVal,

            const Type& correction

        ) -> DivType
        {
            return area & ((lambda*(ownVal - neiVal) + neiVal) + correction);
        };

        fvc::surfaceSum
        (
            this->tinterpScheme_().weights(vf),
            vf,
            faceCorr,
            dotInterpolate,
            div,
            false
        );
    }
    else
    {
        const auto dotInterpolate = [&]
        (
            const vector& area,
            const scalar lambda,
            const Type& ownVal,
            const Type& neiVal
        ) -> DivType
        {
            return area & (lambda*(ownVal - neiVal) + neiVal);
        };

        fvc::surfaceSum
        (
            this->tinterpScheme_().weights(vf),
            vf,
            dotInterpolate,
            div,
            false
        );
    }

    tDiv.ref().primitiveFieldRef() /= mesh.V();

    tDiv.ref().correctBoundaryConditions();

    return tDiv;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

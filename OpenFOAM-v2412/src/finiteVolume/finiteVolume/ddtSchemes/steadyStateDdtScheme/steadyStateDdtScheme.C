/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "steadyStateDdtScheme.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+dt.name()+')',
            mesh().time().timeName(),
            mesh().thisDb()
        ),
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+vf.name()+')',
            mesh().time().timeName(),
            mesh().thisDb()
        ),
        mesh(),
        dimensioned<Type>(vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+rho.name()+','+vf.name()+')',
            mesh().time().timeName(),
            mesh().thisDb()
        ),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+rho.name()+','+vf.name()+')',
            mesh().time().timeName(),
            mesh().thisDb()
        ),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
steadyStateDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
            mesh().time().timeName(),
            mesh().thisDb()
        ),
        mesh(),
        dimensioned<Type>(rho.dimensions()*vf.dimensions()/dimTime, Zero)
    );
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
steadyStateDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );

    return tfvm;
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            mesh(),
            Foam::zero{},  // value
            Uf.dimensions()*dimArea/dimTime
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            mesh(),
            Foam::zero{},  // value
            phi.dimensions()/dimTime
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr("
              + rho.name()
              + ',' + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            mesh(),
            Foam::zero{},  // value
            Uf.dimensions()*dimArea/dimTime
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<typename steadyStateDdtScheme<Type>::fluxFieldType>
steadyStateDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    tmp<fluxFieldType> tCorr
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr("
              + rho.name()
              + ',' + U.name() + ',' + phi.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            mesh(),
            Foam::zero{},  // value
            phi.dimensions()/dimTime
        )
    );

    tCorr.ref().setOriented();

    return tCorr;
}


template<class Type>
tmp<surfaceScalarField> steadyStateDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    auto tphi
    (
        tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "meshPhi",
                mesh().time().timeName(),
                mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedScalar(dimVolume/dimTime, Zero)
        )
    );
    tphi.ref().setOriented();
    return tphi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

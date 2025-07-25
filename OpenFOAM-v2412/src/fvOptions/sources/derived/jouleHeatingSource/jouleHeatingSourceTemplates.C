/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "emptyFvPatchField.H"

template<class Type>
void Foam::fv::jouleHeatingSource::initialiseSigma
(
    const dictionary& dict,
    autoPtr<Function1<Type>>& sigmaVsTPtr
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    IOobject io
    (
        IOobject::scopedName(typeName, "sigma"),
        mesh_.time().timeName(),
        mesh_.thisDb(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        IOobject::REGISTER
    );

    autoPtr<FieldType> tsigma;

    if (dict.found("sigma"))
    {
        // Sigma to be defined using a Function1 type
        sigmaVsTPtr = Function1<Type>::New("sigma", dict, &mesh_);

        tsigma.reset
        (
            new FieldType
            (
                io,
                mesh_,
                Foam::zero{},  // value
                sqr(dimCurrent)/dimPower/dimLength
            )
        );

        Info<< "    Conductivity 'sigma' read from dictionary as f(T)"
            << nl << endl;
    }
    else
    {
        // Sigma to be defined by user input
        io.readOpt(IOobject::MUST_READ);

        tsigma.reset(new FieldType(io, mesh_));

        Info<< "    Conductivity 'sigma' read from file" << nl << endl;
    }

    regIOobject::store(tsigma);
}


template<class Type>
const Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>&
Foam::fv::jouleHeatingSource::updateSigma
(
    const autoPtr<Function1<Type>>& sigmaVsTPtr
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    auto& sigma = mesh_.lookupObjectRef<FieldType>
    (
        IOobject::scopedName(typeName, "sigma")
    );

    if (!sigmaVsTPtr)
    {
        // Electrical conductivity field, sigma, was specified by the user
        return sigma;
    }

    const auto& T = mesh_.lookupObject<volScalarField>(TName_);

    // Internal field
    forAll(sigma, i)
    {
        sigma[i] = sigmaVsTPtr->value(T[i]);
    }


    // Boundary field
    auto& bf = sigma.boundaryFieldRef();
    forAll(bf, patchi)
    {
        fvPatchField<Type>& pf = bf[patchi];
        if (!isA<emptyFvPatchField<Type>>(pf))
        {
            const scalarField& Tbf = T.boundaryField()[patchi];
            forAll(pf, facei)
            {
                pf[facei] = sigmaVsTPtr->value(Tbf[facei]);
            }
        }
    }

    // Update processor patches
    sigma.correctBoundaryConditions();

    return sigma;
}


// ************************************************************************* //

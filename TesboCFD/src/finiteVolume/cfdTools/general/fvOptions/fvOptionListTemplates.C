/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
    Copyright (C) 2020,2023 PCOpt/NTUA
    Copyright (C) 2020,2023 FOSS GP
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

#include "profiling.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::source
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName,
    const dimensionSet& ds
)
{
    checkApplied();

    auto tmtx = tmp<fvMatrix<Type>>::New(field, ds);
    auto& mtx = tmtx.ref();

    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(fvopt, "fvOption().", source.name());

            source.setApplied(fieldi);

            const bool ok = source.isActive();

            if (debug)
            {
                if (ok)
                {
                    Info<< "Apply";
                }
                else
                {
                    Info<< "(Inactive)";
                }
                Info<< " source " << source.name()
                    << " for field " << fieldName << endl;
            }

            if (ok)
            {
                source.addSup(mtx, fieldi);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    return source(field, fieldName, field.dimensions()/dimTime*dimVolume);
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds
    (
        rho.dimensions()*field.dimensions()/dimTime*dimVolume
    );

    auto tmtx = tmp<fvMatrix<Type>>::New(field, ds);
    auto& mtx = tmtx.ref();

    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(fvopt, "fvOption().", source.name());

            source.setApplied(fieldi);

            const bool ok = source.isActive();

            if (debug)
            {
                if (ok)
                {
                    Info<< "Apply";
                }
                else
                {
                    Info<< "(Inactive)";
                }
                Info<< " source " << source.name()
                    << " for field " << fieldName << endl;
            }

            if (ok)
            {
                source.addSup(rho, mtx, fieldi);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(alpha, rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds
    (
        alpha.dimensions()*rho.dimensions()*field.dimensions()
       /dimTime*dimVolume
    );

    auto tmtx = tmp<fvMatrix<Type>>::New(field, ds);
    auto& mtx = tmtx.ref();

    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(fvopt, "fvOption().", source.name());

            source.setApplied(fieldi);

            const bool ok = source.isActive();

            if (debug)
            {
                if (ok)
                {
                    Info<< "Apply";
                }
                else
                {
                    Info<< "(Inactive)";
                }
                Info<< " source " << source.name()
                    << " for field " << fieldName << endl;
            }

            if (ok)
            {
                source.addSup(alpha, rho, mtx, fieldi);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const volScalarField& alpha,
    const geometricOneField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    volScalarField one
    (
        IOobject
        (
            "one",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        this->mesh_,
        dimensionedScalar("one", dimless, scalar(1))
    );

    return this->operator()(alpha, one, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::operator()
(
    const geometricOneField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->operator()(rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::d2dt2
(
    GeometricField<Type, fvPatchField, volMesh>& field
)
{
    return this->d2dt2(field, field.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionList::d2dt2
(
    GeometricField<Type, fvPatchField, volMesh>& field,
    const word& fieldName
)
{
    return source(field, fieldName, field.dimensions()/sqr(dimTime)*dimVolume);
}


template<class Type>
void Foam::fv::optionList::constrain(fvMatrix<Type>& eqn)
{
    checkApplied();

    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(eqn.psi().name());

        if (fieldi != -1)
        {
            addProfiling(fvopt, "fvOption::constrain.", eqn.psi().name());

            source.setApplied(fieldi);

            const bool ok = source.isActive();

            if (debug)
            {
                if (ok)
                {
                    Info<< "Constrain";
                }
                else
                {
                    Info<< "(Inactive constrain)";
                }
                Info<< " source " << source.name()
                    << " for field " << eqn.psi().name() << endl;
            }

            if (ok)
            {
                source.constrain(eqn, fieldi);
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fv::optionList::correct
(
    GeometricField<Type, PatchField, GeoMesh>& field
)
{
    const word& fieldName = field.name();

    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(fvopt, "fvOption::correct.", source.name());

            source.setApplied(fieldi);

            const bool ok = source.isActive();

            if (debug)
            {
                if (ok)
                {
                    Info<< "Correct";
                }
                else
                {
                    Info<< "(Inactive correct)";
                }
                Info<< " source " << source.name()
                    << " for field " << fieldName << endl;
            }

            if (ok)
            {
                source.correct(field);
            }
        }
    }
}


template<class Type>
void Foam::fv::optionList::postProcessSens
(
    Field<Type>& sensField,
    const word& fieldName,
    const word& designVariablesName
)
{
    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(fvopt, "fvOption::postProcessSens.", source.name());

            const bool ok = source.isActive();

            if (debug && ok)
            {
                Info<< "Post processing sensitivity source "
                    << source.name() << " for field " << fieldName << endl;
            }

            if (ok)
            {
                source.postProcessSens
                (
                    sensField,
                    fieldName,
                    designVariablesName
                );
            }
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fv::optionList::postProcessAuxSens
(
    const GeometricField<Type, PatchField, GeoMesh>& primal,
    const GeometricField<Type, PatchField, GeoMesh>& adjoint,
    scalarField& sensField,
    const word& fieldName
)
{
    for (fv::option& source : *this)
    {
        const label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling
            (
                fvopt,
                "fvOption::postProcessAuxSens.", source.name()
            );

            const bool ok = source.isActive();

            if (debug && ok)
            {
                Info<< "Post processing sensitivity source "
                    << source.name() << " for field " << fieldName << endl;
            }

            if (ok)
            {
                source.postProcessAuxSens
                (
                    primal,
                    adjoint,
                    sensField,
                    fieldName
                );
            }
        }
    }
}


// ************************************************************************* //

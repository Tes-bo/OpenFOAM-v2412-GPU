/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SRFModel.H"
#include "SRFVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace SRF
    {
        defineTypeNameAndDebug(SRFModel, 0);
        defineRunTimeSelectionTable(SRFModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SRF::SRFModel::SRFModel
(
    const word& type,
    const volVectorField& Urel
)
:
    IOdictionary
    (
        IOobject
        (
            "SRFProperties",
            Urel.time().constant(),
            Urel.db(),
            IOobject::READ_MODIFIED,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        )
    ),
    Urel_(Urel),
    mesh_(Urel_.mesh()),
    origin_("origin", dimLength, get<vector>("origin")),
    axis_(normalised(get<vector>("axis"))),
    SRFModelCoeffs_(optionalSubDict(type + "Coeffs")),
    omega_(dimensionedVector("omega", dimless/dimTime, Zero))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SRF::SRFModel::~SRFModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SRF::SRFModel::read()
{
    if (regIOobject::read())
    {
        // Re-read origin
        readEntry("origin", origin_);

        // Re-read axis
        readEntry("axis", axis_);
        axis_.normalise();

        // Re-read sub-model coeffs
        SRFModelCoeffs_ = optionalSubDict(type() + "Coeffs");

        return true;
    }

    return false;
}


const Foam::dimensionedVector& Foam::SRF::SRFModel::origin() const
{
    return origin_;
}


const Foam::vector& Foam::SRF::SRFModel::axis() const
{
    return axis_;
}


const Foam::dimensionedVector& Foam::SRF::SRFModel::omega() const
{
    return omega_;
}


Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::SRF::SRFModel::Fcoriolis() const
{
    return volVectorField::Internal::New
    (
        "Fcoriolis",
        IOobject::NO_REGISTER,
        (
            2.0*omega_ ^ Urel_.internalField()
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::SRF::SRFModel::Fcentrifugal() const
{
    return volVectorField::Internal::New
    (
        "Fcentrifugal",
        IOobject::NO_REGISTER,
        (
            omega_ ^ (omega_ ^ (mesh_.C().internalField() - origin_))
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::vector, Foam::volMesh>>
Foam::SRF::SRFModel::Su() const
{
    return Fcoriolis() + Fcentrifugal();
}


Foam::vectorField Foam::SRF::SRFModel::velocity
(
    const vectorField& positions
) const
{
    return vectorField
    (
        omega_.value()
      ^ (
            (positions - origin_.value())
          - axis_*(axis_ & (positions - origin_.value()))
        )
    );
}


Foam::tmp<Foam::volVectorField> Foam::SRF::SRFModel::U() const
{
    const int oldLocal = volVectorField::Boundary::localConsistency;
    volVectorField::Boundary::localConsistency = 0;
    tmp<volVectorField> relPos(mesh_.C() - origin_);

    auto tU = volVectorField::New
    (
        "Usrf",
        IOobject::NO_REGISTER,
        (
            omega_ ^ (relPos() - axis_*(axis_ & relPos()))
        )
    );
    volVectorField::Boundary::localConsistency = oldLocal;

    return tU;
}


Foam::tmp<Foam::volVectorField> Foam::SRF::SRFModel::Uabs() const
{
    tmp<volVectorField> Usrf = U();

    auto tUabs = volVectorField::New
    (
        "Uabs",
        IOobject::NO_REGISTER,
        Usrf
    );
    auto& Uabs = tUabs.ref();

    // Add SRF contribution to internal field
    Uabs.primitiveFieldRef() += Urel_.primitiveField();

    // Add Urel boundary contributions
    volVectorField::Boundary& Uabsbf = Uabs.boundaryFieldRef();
    const volVectorField::Boundary& bvf = Urel_.boundaryField();

    forAll(bvf, i)
    {
        if (isA<SRFVelocityFvPatchVectorField>(bvf[i]))
        {
            // Only include relative contributions from
            // SRFVelocityFvPatchVectorField's
            const SRFVelocityFvPatchVectorField& UrelPatch =
                refCast<const SRFVelocityFvPatchVectorField>(bvf[i]);
            if (UrelPatch.relative())
            {
                Uabsbf[i] += Urel_.boundaryField()[i];
            }
        }
        else
        {
            Uabsbf[i] += Urel_.boundaryField()[i];
        }
    }

    return tUabs;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "sampledSurface.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurface, 0);
    defineRunTimeSelectionTable(sampledSurface, word);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::sampledSurface::clearGeom() const
{
    area_ = -1;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sampledSurface> Foam::sampledSurface::New
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const word sampleType(dict.get<word>("type"));

    DebugInfo
        << "Selecting sampledType " << sampleType << endl;

    auto* ctorPtr = wordConstructorTable(sampleType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "sample",
            sampleType,
            *wordConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<sampledSurface>(ctorPtr(name, mesh, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurface::sampledSurface(const word& name, std::nullptr_t)
:
    name_(name),
    mesh_(NullObjectRef<polyMesh>()),
    enabled_(true),
    invariant_(false),
    isPointData_(false),
    area_(-1)
{}


Foam::sampledSurface::sampledSurface
(
    const word& name,
    const polyMesh& mesh,
    const bool interpolateToPoints
)
:
    name_(name),
    mesh_(mesh),
    enabled_(true),
    invariant_(false),
    isPointData_(interpolateToPoints),
    area_(-1)
{}


Foam::sampledSurface::sampledSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(dict.getOrDefault<word>("name", name)),
    mesh_(mesh),
    enabled_(dict.getOrDefault("enabled", true)),
    invariant_(dict.getOrDefault("invariant", false)),
    isPointData_(dict.getOrDefault("interpolate", false)),
    area_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurface::~sampledSurface()
{
    clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::sampledSurface::area() const
{
    if (area_ < 0)
    {
        area_ = gSum(magSf());
    }

    return area_;
}


bool Foam::sampledSurface::isPointData(const bool on)
{
    bool old(isPointData_);
    isPointData_ = on;
    return old;
}


bool Foam::sampledSurface::withSurfaceFields() const
{
    return false;
}


Foam::tmp<Foam::scalarField> Foam::sampledSurface::sample
(
    const surfaceScalarField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::vectorField> Foam::sampledSurface::sample
(
    const surfaceVectorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurface::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurface::sample
(
    const surfaceSymmTensorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::tensorField> Foam::sampledSurface::sample
(
    const surfaceTensorField& sField
) const
{
    NotImplemented;
    return nullptr;
}


void Foam::sampledSurface::print(Ostream& os, int level) const
{
    os << type();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const sampledSurface& s)
{
    // Print with more information
    s.print(os, 1);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

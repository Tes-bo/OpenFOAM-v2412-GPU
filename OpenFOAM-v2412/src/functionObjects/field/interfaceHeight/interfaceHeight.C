/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "interfaceHeight.H"
#include "fvMesh.H"
#include "interpolation.H"
#include "IOmanip.H"
#include "meshSearch.H"
#include "midPointAndFaceSet.H"
#include "Time.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(interfaceHeight, 0);
    addToRunTimeSelectionTable(functionObject, interfaceHeight, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceHeight::writePositions()
{
    const uniformDimensionedVectorField& g =
        mesh_.time().lookupObject<uniformDimensionedVectorField>("g");
    vector gHat = vector::zero;

    if (mag(direction_) > 0.0)
    {
        gHat = direction_/mag(direction_);
    }
    else
    {
        gHat = g.value()/mag(g.value());
    }

    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    autoPtr<interpolation<scalar>>
        interpolator
        (
            interpolation<scalar>::New(interpolationScheme_, alpha)
        );

    if (Pstream::master())
    {
        files(fileID::heightFile) << mesh_.time().timeName() << tab;
        files(fileID::positionFile) << mesh_.time().timeName() << tab;
    }

    forAll(locations_, li)
    {
        // Create a set along a ray projected in the direction of gravity
        const midPointAndFaceSet set
        (
            "",
            mesh_,
            meshSearch(mesh_),
            "xyz",
            locations_[li] + gHat*mesh_.bounds().mag(),
            locations_[li] - gHat*mesh_.bounds().mag()
        );

        // Find the height of the location above the boundary
        scalar hLB = set.size() ? - gHat & (locations_[li] - set[0]) : - GREAT;
        reduce(hLB, maxOp<scalar>());

        // Calculate the integrals of length and length*alpha along the sampling
        // line. The latter is equal to the equivalent length with alpha equal
        // to one.
        scalar sumLength = 0, sumLengthAlpha = 0;
        for(label si = 0; si < set.size() - 1; ++ si)
        {
            if (set.segments()[si] != set.segments()[si+1])
            {
                continue;
            }

            const vector& p0 = set[si], p1 = set[si+1];
            const label c0 = set.cells()[si], c1 = set.cells()[si+1];
            const label f0 = set.faces()[si], f1 = set.faces()[si+1];
            const scalar a0 = interpolator->interpolate(p0, c0, f0);
            const scalar a1 = interpolator->interpolate(p1, c1, f1);

            const scalar l = - gHat & (p1 - p0);
            sumLength += l;
            sumLengthAlpha += l*(a0 + a1)/2;
        }

        reduce(sumLength, sumOp<scalar>());
        reduce(sumLengthAlpha, sumOp<scalar>());

        // Write out
        if (Pstream::master())
        {
            // Interface heights above the boundary and location
            const scalar hIB =
                liquid_ ? sumLengthAlpha : sumLength - sumLengthAlpha;
            const scalar hIL = hIB - hLB;

            // Position of the interface
            const point p = locations_[li] - gHat*hIL;

            const Foam::Omanip<int> w = valueWidth(1);

            files(fileID::heightFile) << w << hIB << w << hIL;
            files(fileID::positionFile) << '(' << w << p.x() << w << p.y()
                << valueWidth() << p.z() << ") ";
        }
    }

    if (Pstream::master())
    {
        files(fileID::heightFile).endl();
        files(fileID::positionFile).endl();
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::interfaceHeight::writeFileHeader(const fileID fid)
{
    forAll(locations_, li)
    {
        writeHeaderValue
        (
            files(fid),
            "Location " + Foam::name(li),
            locations_[li]
        );
    }

    switch (fileID(fid))
    {
        case fileID::heightFile:
        {
            writeHeaderValue
            (
                files(fid),
                "hB",
                "Interface height above the boundary"
            );
            writeHeaderValue
            (
                files(fid),
                "hL",
                "Interface height above the location"
            );
            break;
        }
        case fileID::positionFile:
        {
            writeHeaderValue(files(fid), "p", "Interface position");
            break;
        }
    }

    const Foam::Omanip<int> w = valueWidth(1);

    writeCommented(files(fid), "Location");
    forAll(locations_, li)
    {
        switch (fid)
        {
            case fileID::heightFile:
                files(fid) << w << li << w << ' ';
                break;
            case fileID::positionFile:
                files(fid) << w << li << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    files(fid).endl();

    writeCommented(files(fid), "Time");
    forAll(locations_, li)
    {
        switch (fid)
        {
            case fileID::heightFile:
                files(fid) << w << "hB" << w << "hL";
                break;
            case fileID::positionFile:
                files(fid) << w << "p" << w << ' ' << w << ' ' << "  ";
                break;
        }
    }
    files(fid).endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::interfaceHeight::interfaceHeight
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    liquid_(true),
    alphaName_("alpha"),
    interpolationScheme_("cellPoint"),
    direction_(vector::zero),
    locations_()
{
    read(dict);
    resetNames({"height", "position"});

    if (Pstream::master())
    {
        writeFileHeader(fileID::heightFile);
        writeFileHeader(fileID::positionFile);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::interfaceHeight::read(const dictionary& dict)
{
    dict.readIfPresent("alpha", alphaName_);
    dict.readIfPresent("liquid", liquid_);
    dict.readEntry("locations", locations_);
    dict.readIfPresent("interpolationScheme", interpolationScheme_);
    dict.readIfPresent("direction", direction_);

    return true;
}


bool Foam::functionObjects::interfaceHeight::execute()
{
    return true;
}


bool Foam::functionObjects::interfaceHeight::end()
{
    return true;
}


bool Foam::functionObjects::interfaceHeight::write()
{
    logFiles::write();

    writePositions();

    return true;
}


// ************************************************************************* //

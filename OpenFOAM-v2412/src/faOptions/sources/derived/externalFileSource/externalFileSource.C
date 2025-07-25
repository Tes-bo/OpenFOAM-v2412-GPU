/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "externalFileSource.H"
#include "fam.H"
#include "faScalarMatrix.H"
#include "zeroGradientFaPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(externalFileSource, 0);
    addToRunTimeSelectionTable(option, externalFileSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::externalFileSource::externalFileSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& m
)
:
    fa::faceSetOption(sourceName, modelType, dict, m),
    fieldName_(dict.get<word>("fieldName")),
    tableName_(dict.get<word>("tableName")),
    pExt_
    (
        IOobject
        (
            "pExt",
            regionMesh().thisDb().time().timeName(),
            regionMesh().thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            (
                dict.getOrDefault("store", false)
              ? IOobject::REGISTER
              : IOobject::NO_REGISTER
            )
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    curTimeIndex_(-1),
    mapping_()
{
    fieldNames_.resize(1, fieldName_);

    /// FUTURE?
    /// // Optional entry (mandatory = false)
    /// pExt_.dimensions().readEntry("dimensions",  dict, false);

    fa::option::resetApplied();

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::externalFileSource::updateMapping()
{
    // Set up mapped values per patch
    const scalar t = mesh().time().value();

    PtrList<scalarField> patchValues(mapping_.size());

    forAll(mapping_, patchi)
    {
        const auto* map = mapping_.get(patchi);

        if (map)
        {
            patchValues.set(patchi, map->value(t));
        }
    }

    vsm().mapToSurface<scalar>(patchValues, pExt_.field());

    // Zero pressure for non-mapped faces
    faceSetOption::subsetFilter(pExt_.field());
}


void Foam::fa::externalFileSource::addSup
(
    const areaScalarField& solidMass,
    faMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isActive())
    {
        DebugInfo
            << name() << ": applying source to "
            << eqn.psi().name() << endl;

        if (curTimeIndex_ != mesh().time().timeIndex())
        {
            updateMapping();

            eqn += pExt_/solidMass.internalField();

            curTimeIndex_ = mesh().time().timeIndex();
        }
    }
}


bool Foam::fa::externalFileSource::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        // Set up mapping (per-patch) for referenced polyPatches (sorted order)
        // - size is maxPolyPatch+1

        const labelList& patches = regionMesh().whichPolyPatches();

        mapping_.clear();
        mapping_.resize(patches.empty() ? 0 : (patches.back()+1));

        for (const label patchi : patches)
        {
            const polyPatch& p = mesh_.boundaryMesh()[patchi];

            mapping_.set
            (
                patchi,
                (
                    new PatchFunction1Types::MappedFile<scalar>
                    (
                        p,
                        "uniformValue", // entryName
                        dict,
                        tableName_,     // field table name
                        true            // face values
                    )
                )
            );
        }

        return true;
    }

    return false;
}


// ************************************************************************* //

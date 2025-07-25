/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "patchCellsSource.H"
#include "boundarySourcePatch.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(patchCellsSource, 0);
    addToRunTimeSelectionTable(option, patchCellsSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::patchCellsSource::patchCellsSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(sourceName, modelType, dict, mesh),
    curTimeIndex_(-1),
    isEnergySource_(false)
{
    fieldNames_.resize(1);

    label nFields = 0;

    if
    (
        coeffs_.readIfPresent("U", fieldNames_[0])
     && fieldNames_[0] != "none"
    )
    {
        ++nFields;
    }

    if
    (
        coeffs_.readIfPresent("he", fieldNames_[0])
     && fieldNames_[0] != "none"
    )
    {
        isEnergySource_ = true;
        ++nFields;
    }

    if
    (
        coeffs_.readIfPresent("species", fieldNames_[0])
     && fieldNames_[0] != "none"
    )
    {
        ++nFields;
    }

    if (nFields != 1)
    {
        FatalIOErrorInFunction(coeffs_)
            << "Must be specified for one field (U | he | species), but "
            << nFields << " fields were specified!" << endl
            << exit(FatalIOError);
    }

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::patchCellsSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (curTimeIndex_ == mesh_.time().timeIndex())
    {
        return;
    }
    curTimeIndex_ = mesh_.time().timeIndex();

    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    auto tsu = DimensionedField<scalar, volMesh>::New
    (
        name_ + eqn.psi().name() + "Su",
        mesh_,
        dimensionedScalar(eqn.dimensions()/dimVolume, Zero)
    );

    // If source applied to he, we need to loop over T for BC's
    const volScalarField& psi =
    (
        isEnergySource_
      ? mesh_.lookupObject<volScalarField>("T")
      : mesh_.lookupObject<volScalarField>(eqn.psi().name())
    );
    const auto& psib = psi.boundaryField();

    forAll(psib, patchi)
    {
        const auto* bpatchPtr = isA<boundarySourcePatch>(psib[patchi]);

        if (bpatchPtr)
        {
            tmp<scalarField> tsbnd = bpatchPtr->patchSource();
            const auto& sbnd = tsbnd.cref();

            UIndirectList<scalar>
            (
                tsu.ref().field(),
                mesh_.boundary()[patchi].faceCells()
            ) = sbnd;
        }
    }

    if (debug)
    {
         Info<< "Field source rate min/max : " << gMinMax(tsu()) << endl;
    }

    eqn += tsu;
}


bool Foam::fv::patchCellsSource::read(const dictionary& dict)
{
    if (!fv::option::read(dict))
    {
        return false;
    }

    return true;
}


// ************************************************************************* //

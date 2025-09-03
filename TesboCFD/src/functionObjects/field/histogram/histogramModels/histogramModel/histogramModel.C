/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "histogramModel.H"
#include "fvMesh.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(histogramModel, 0);
    defineRunTimeSelectionTable(histogramModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::histogramModel::writeFileHeader(Ostream& os)
{
    writeHeader(os, "Histogram");
    writeCommented(os, "Time");
    writeTabbed(os, "binMidPoints");
    writeTabbed(os, "dataCounts");
    writeTabbed(os, "dataValues");
    os  << endl;
}


Foam::volScalarField& Foam::histogramModel::getOrReadField
(
    const word& fieldName
) const
{
    auto* ptr = mesh_.getObjectPtr<volScalarField>(fieldName);

    if (!ptr)
    {
        ptr = new volScalarField
        (
            IOobject
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_.thisDb(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                IOobject::REGISTER
            ),
            mesh_
        );
        regIOobject::store(ptr);
    }

    return *ptr;
}


void Foam::histogramModel::write
(
    scalarField& dataNormalised,
    const labelField& dataCount,
    const scalarField& magBinMidPoint
)
{
    if (!Pstream::master())
    {
        return;
    }

    const scalar sumData = sum(dataNormalised);

    if (sumData < SMALL)
    {
        return;
    }

    dataNormalised /= sumData;

    const auto time = mesh().time().value();

    forAll(dataNormalised, i)
    {
        file()
            << time << tab
            << magBinMidPoint[i] << tab
            << dataCount[i] << tab
            << dataNormalised[i]
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::histogramModel::histogramModel
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    writeFile(mesh, name, "histogram", dict),
    mesh_(mesh),
    fieldName_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::histogramModel::read(const dictionary& dict)
{
    if (!functionObjects::writeFile::read(dict))
    {
        return false;
    }

    fieldName_ = dict.get<word>("field");

    if (writeToFile() && !writtenHeader_)
    {
        writeFileHeader(file());
    }

    return true;
}


// ************************************************************************* //

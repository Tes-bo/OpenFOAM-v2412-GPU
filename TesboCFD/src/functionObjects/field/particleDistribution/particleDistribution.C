/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "particleDistribution.H"
#include "addToRunTimeSelectionTable.H"
#include "general.H"
#include "fvMesh.H"
#include "cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(particleDistribution, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        particleDistribution,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::particleDistribution::particleDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(runTime, name),
    cloudName_("unknown-cloudName"),
    tagFieldName_("none"),
    rndGen_(),
    nameVsBinWidth_(),
    writerPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::particleDistribution::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.readEntry("cloud", cloudName_);
        dict.readIfPresent("tagField", tagFieldName_);
        dict.readEntry("nameVsBinWidth", nameVsBinWidth_);

        const word setFormat(dict.get<word>("setFormat"));
        writerPtr_ = coordSetWriter::New
        (
            setFormat,
            dict.subOrEmptyDict("formatOptions").optionalSubDict(setFormat)
        );

        Info<< type() << " " << name() << " output:" << nl
            << "    Processing cloud : " << cloudName_ << nl
            << endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::particleDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::particleDistribution::write()
{
    Log << type() << " " << name() << " output:" << endl;

    const cloud* cloudPtr = mesh_.cfindObject<cloud>(cloudName_);

    if (!cloudPtr)
    {
        WarningInFunction
            << "Unable to find cloud " << cloudName_
            << " in the mesh database.  Available clouds include:"
            << flatOutput(mesh_.sortedNames<cloud>()) << endl;

        return false;
    }

    const cloud& c = *cloudPtr;

    objectRegistry cloudObr
    (
        IOobject
        (
            scopedName("CloudRegistry"),
            mesh_.time().timeName(),
            cloud::prefix,
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    c.writeObjects(cloudObr);

    List<DynamicList<label>> tagAddr;
    if
    (
        tagFieldName_ != "none"
     && cloudObr.foundObject<IOField<scalar>>(tagFieldName_)
    )
    {
        // Tag field present - generate distribution per tag
        const IOField<label>& tag =
            cloudObr.lookupObject<IOField<label>>(tagFieldName_);
        const labelHashSet tagMap(tag);
        const label tagMax = tagMap.size();

        List<DynamicList<label>> tagAddr(tagMax);
        forAll(tag, i)
        {
            label newTag = tagMap[tag[i]];
            tagAddr[newTag].append(i);
        }
    }


    forAll(nameVsBinWidth_, i)
    {
        const bool ok
        (
            processField<scalar>(cloudObr, i, tagAddr)
         || processField<vector>(cloudObr, i, tagAddr)
         || processField<tensor>(cloudObr, i, tagAddr)
         || processField<sphericalTensor>(cloudObr, i, tagAddr)
         || processField<symmTensor>(cloudObr, i, tagAddr)
         || processField<tensor>(cloudObr, i, tagAddr)
        );

        if (log && !ok)
        {
            WarningInFunction
                << "Unable to find field " << nameVsBinWidth_[i].first()
                << " in the " << cloudName_ << " cloud database" << endl;
        }
    }

    return true;
}


void Foam::functionObjects::particleDistribution::generateDistribution
(
    const word& fieldName,
    const scalarField& field,
    const scalar binWidth,
    const label tag
)
{
    if (field.empty())
    {
        return;
    }

    word fldName(fieldName);
    if (tag != -1)
    {
        fldName += '_' + Foam::name(tag);
    }

    distributionModels::general distribution
    (
        field,
        binWidth,
        rndGen_
    );

    Field<scalar> distX(distribution.x());
    Field<scalar> distY(distribution.y());

    pointField xBin(distX.size(), Zero);
    xBin.replace(vector::X, distX);

    const coordSet coords(fldName, "x", std::move(xBin), std::move(distX));

    writerPtr_->open(coords, baseTimeDir() / fldName);
    fileName outFile = writerPtr_->write(fldName, distY);
    writerPtr_->close(true);

    Log << "    Wrote distribution of " << fieldName
        << " to " << time_.relativePath(outFile) << endl;
}


// ************************************************************************* //

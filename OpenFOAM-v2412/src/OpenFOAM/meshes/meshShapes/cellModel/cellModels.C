/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "cellModel.H"
#include "etcFiles.H"
#include "IFstream.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::PtrList<Foam::cellModel> Foam::cellModel::models_;

Foam::List<const Foam::cellModel*> Foam::cellModel::modelPtrs_;

const Foam::Enum<Foam::cellModel::modelType> Foam::cellModel::modelNames
({
    { modelType::UNKNOWN, "unknown" },
    { modelType::HEX, "hex" },
    { modelType::WEDGE, "wedge" },
    { modelType::PRISM, "prism" },
    { modelType::PYR, "pyr" },
    { modelType::TET, "tet" },
    { modelType::TETWEDGE, "tetWedge" },
    { modelType::SPLITHEX, "splitHex" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::cellModel::constructModels()
{
    if (models_.size())
    {
        FatalErrorInFunction
            << "attempt to re-construct cellModeller when it already exists"
            << exit(FatalError);
    }

    IFstream is(findEtcFile("cellModels", true));  // Mandatory file

    PtrList<cellModel> newPtrs(is);
    models_.swap(newPtrs);

    ///Info<< "loading " << models_.size()
    ///    << " cell models from etc/controlDict" << endl;


    // Build two lookups: by index, by name
    // Since there are relatively few models, use straight lookup for the index
    // and a linear (non-hashed) search for the name.
    // Lookup by name is less likely than lookup by enum anyhow.

    label maxIndex = 0;
    forAll(models_, i)
    {
        if (maxIndex < models_[i].index())
        {
            maxIndex = models_[i].index();
        }
    }

    modelPtrs_.clear();
    modelPtrs_.setSize(maxIndex+1, nullptr);

    wordHashSet used(2*maxIndex);

    forAll(models_, i)
    {
        const label modelIndex = models_[i].index();
        const word& modelName = models_[i].name();
        const cellModel* ptr = &models_[i];

        if (used.insert(modelName))
        {
            if (modelPtrs_[modelIndex])
            {
                FatalErrorInFunction
                    << "more than one model share the index "
                    << modelIndex
                    << exit(FatalError);
            }

            modelPtrs_[modelIndex] = ptr;
        }
        else
        {
            FatalErrorInFunction
                << "more than one model share the name "
                << modelName
                << exit(FatalError);
        }
    }
}


const Foam::cellModel* Foam::cellModel::ptr(const modelType model)
{
    return ptr(label(model));
}


const Foam::cellModel* Foam::cellModel::ptr(const word& modelName)
{
    if (models_.empty())
    {
        constructModels();
    }

    const label n = models_.size();
    for (label i = 0; i < n; ++i)
    {
        if (models_[i].name() == modelName)
        {
            return &(models_[i]);
        }
    }

    return nullptr;
}


const Foam::cellModel* Foam::cellModel::ptr(const label modelIndex)
{
    if (models_.empty())
    {
        constructModels();
    }

    return (modelIndex < modelPtrs_.size() ? modelPtrs_[modelIndex] : nullptr);
}


const Foam::cellModel& Foam::cellModel::ref(const modelType model)
{
    const cellModel* p = ptr(model);

    if (!p)
    {
        FatalErrorInFunction
            << "No such cellModel: " << modelNames[model]
            << exit(FatalError);
    }

    return *p;
}


const Foam::cellModel& Foam::cellModel::ref(const word& modelName)
{
    const cellModel* p = ptr(modelName);

    if (!p)
    {
        FatalErrorInFunction
            << "No such cellModel: " << modelName
            << exit(FatalError);
    }

    return *p;
}


const Foam::cellModel& Foam::cellModel::ref(const label modelIndex)
{
    const cellModel* p = ptr(modelIndex);

    if (!p)
    {
        FatalErrorInFunction
            << "No such cellModel: " << modelIndex
            << exit(FatalError);
    }

    return *p;
}


// ************************************************************************* //

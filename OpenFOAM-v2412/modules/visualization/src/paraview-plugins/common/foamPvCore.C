/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "foamPvCore.H"

#include "memInfo.H"
#include "DynamicList.H"

#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkInformation.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::Ostream& Foam::foamPvCore::printDataArraySelection
(
    Ostream& os,
    vtkDataArraySelection* select
)
{
    if (!select)
    {
        return os;
    }

    const int n = select->GetNumberOfArrays();

    os << n << '(';
    for (int i=0; i < n; ++i)
    {
        if (i) os << ' ';
        os  << select->GetArrayName(i) << '='
            << (select->GetArraySetting(i) ? 1 : 0);
    }
    os << ')';

    return os;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamPvCore::addToBlock
(
    vtkMultiBlockDataSet* output,
    vtkDataSet* dataset,
    const arrayRange& selector,
    const label datasetNo,
    const std::string& datasetName
)
{
    const int blockNo = selector.block();

    vtkDataObject* dataObj = output->GetBlock(blockNo);
    vtkSmartPointer<vtkMultiBlockDataSet> block =
        vtkMultiBlockDataSet::SafeDownCast(dataObj);

    if (!block)
    {
        if (dataObj)
        {
            FatalErrorInFunction
                << "Block already has a vtkDataSet assigned to it"
                << endl;
            return;
        }

        block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        output->SetBlock(blockNo, block);
    }

    #ifdef FULLDEBUG
    {
        Info<< "block[" << blockNo << "] has "
            << block->GetNumberOfBlocks()
            <<  " datasets prior to adding set " << datasetNo
            <<  " with name: " << datasetName << endl;
    }
    #endif

    block->SetBlock(datasetNo, dataset);

    // Name the output block when assigning dataset 0
    if (datasetNo == 0)
    {
        output->GetMetaData(blockNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            selector.name()
        );
    }

    if (datasetName.size())
    {
        block->GetMetaData(datasetNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            datasetName.c_str()
        );
    }
}


Foam::wordHashSet Foam::foamPvCore::getSelected
(
    vtkDataArraySelection* select,
    const bool on
)
{
    const int n = select->GetNumberOfArrays();
    wordHashSet selected(2*n);

    for (int i=0; i < n; ++i)
    {
        if (on ? select->GetArraySetting(i) : true)
        {
            selected.set(getFoamName(select->GetArrayName(i)));
        }
    }

    return selected;
}


Foam::wordHashSet Foam::foamPvCore::getSelected
(
    vtkDataArraySelection* select,
    const labelRange& selector,
    const bool on
)
{
    const int n = select->GetNumberOfArrays();
    wordHashSet selected(2*n);

    for (auto i : selector)
    {
        if (on ? select->GetArraySetting(i) : true)
        {
            selected.set(getFoamName(select->GetArrayName(i)));
        }
    }

    return selected;
}


Foam::HashSet<Foam::string>
Foam::foamPvCore::getSelectedArraySet
(
    vtkDataArraySelection* select,
    const bool on
)
{
    const int n = select->GetNumberOfArrays();
    HashSet<string> selected(2*n);

    for (int i=0; i < n; ++i)
    {
        if (on ? select->GetArraySetting(i) : true)
        {
            selected.insert(select->GetArrayName(i));
        }
    }

    #ifdef FULLDEBUG
    if (select)
    {
        Info<< "available(";
        for (int i=0; i < n; ++i)
        {
            Info<< " \"" << select->GetArrayName(i) << "\"";
        }
        Info<< " )\nselected(";

        for (auto k : selected)
        {
            Info<< " " << k;
        }
        Info<< " )\n";
    }
    #endif

    return selected;
}


Foam::Map<Foam::string>
Foam::foamPvCore::getSelectedArrayMap
(
    vtkDataArraySelection* select,
    const bool on
)
{
    const int n = select->GetNumberOfArrays();
    Map<string> selected(2*n);

    for (int i=0; i < n; ++i)
    {
        if (on ? select->GetArraySetting(i) : true)
        {
            selected.insert(i, select->GetArrayName(i));
        }
    }

    return selected;
}


Foam::word Foam::foamPvCore::getFoamName(const std::string& str)
{
    if (str.size())
    {
        auto beg = str.rfind('/');
        if (beg == std::string::npos)
        {
            beg = 0;
        }
        else
        {
            ++beg;
        }

        auto end = beg;

        while (str[end] && word::valid(str[end]))
        {
            ++end;
        }

        // Already checked for valid/invalid chars
        return word(str.substr(beg, beg+end), false);
    }

    return word::null;
}


void Foam::foamPvCore::printMemory()
{
    Info<< "mem peak/size/rss: " << Foam::memInfo{} << endl;
}


// ************************************************************************* //

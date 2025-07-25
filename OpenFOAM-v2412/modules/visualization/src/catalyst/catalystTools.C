/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "catalystTools.H"
#include "Time.H"

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkDataObject.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::wordHashSet Foam::catalyst::dataQuery::emptyWordHash_;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::catalyst::printInfo(Ostream& os, vtkCPDataDescription* descrip)
{
    if (!descrip)
    {
        return;
    }

    const unsigned nItems = descrip->GetNumberOfInputDescriptions();

    for (unsigned itemi = 0; itemi < nItems; ++itemi)
    {
        vtkCPInputDataDescription* input = descrip->GetInputDescription(itemi);
        if (!input) continue;  // should not happen

        os <<"input: " << descrip->GetInputDescriptionName(itemi) << nl;

        const unsigned nFields = input->GetNumberOfFields();
        for (unsigned fieldi = 0; fieldi < nFields; ++fieldi)
        {
            os << "    field: " << input->GetFieldName(fieldi) << nl;
        }
        if (!nFields) os <<"     no fields requested" << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyst::timeQuery::timeQuery(const Foam::Time& when)
:
    timeValue(when.timeOutputValue()),
    timeIndex(when.timeIndex()),
    forced(when.timeOutputValue() >= when.endTime().value())
{}


Foam::catalyst::dataQuery::dataQuery(const Foam::Time& when)
:
    timeQuery(when),
    names_(),
    fields_(),
    active_(),
    descrip_(vtkCPDataDescription::New())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::catalyst::dataQuery::~dataQuery()
{
    descrip_->Delete();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::catalyst::dataQuery::set
(
    const word& channel,
    const wordHashSet& fields
)
{
    auto iter = fields_.find(channel);

    if (iter.good())
    {
        *iter = fields;
    }
    else
    {
        active_.insert(names_.size());
        fields_.insert(channel, fields);
        names_.append(channel);
    }
}


Foam::label Foam::catalyst::dataQuery::query(vtkCPProcessor* coproc)
{
    dataQuery& dataq = *this;

    if (!coproc)
    {
        Info<< "No ParaView Catalyst initialized" << endl;
        dataq.clear();
        return dataq.size();
    }

    if (dataq.empty())
    {
        // No channels names have been published by the simulation
        return dataq.size();
    }

    vtkCPDataDescription* descrip = dataq.description();

    descrip->SetTimeData(dataq.timeValue, dataq.timeIndex);
    descrip->SetForceOutput(dataq.forced);

    const wordList inputNames(dataq.names());

    // Sort out which channels already exist, are new, or disappeared

    // The currently defined channels
    wordHashSet currChannels;

    const unsigned n = descrip->GetNumberOfInputDescriptions();
    for (unsigned i=0; i < n; ++i)
    {
        currChannels.insert
        (
            word::validate(descrip->GetInputDescriptionName(i))
        );
    }

    wordHashSet oldChannels(currChannels);
    oldChannels.erase(inputNames);

    if (oldChannels.size())
    {
        // Some channels disappear - remove and redo everything
        descrip->ResetAll();
        currChannels.clear();
    }

    // Add channels.

    // Note: this misses updating field information for previously
    // existing inputs.

    for (const word& channel : inputNames)
    {
        if (currChannels.contains(channel))
        {
            continue;
        }

        descrip->AddInput(channel.c_str());
        auto* input = descrip->GetInputDescriptionByName(channel.c_str());

        for (const word& fldName : dataq.fields(channel))
        {
            input->AddField(fldName.c_str(), vtkDataObject::CELL);
            input->AddField(fldName.c_str(), vtkDataObject::POINT);
        }
    }

    if
    (
        !coproc->RequestDataDescription(descrip)
     || !descrip->GetIfAnyGridNecessary()
    )
    {
        dataq.clear();
        return dataq.size();
    }

    for (const word& channel : inputNames)
    {
        auto* input = descrip->GetInputDescriptionByName(channel.c_str());

        if (input && input->GetIfGridIsNecessary())
        {
            wordHashSet requestedFields;

            for (const word& fldName : dataq.fields(channel))
            {
                if
                (
                    input->IsFieldNeeded(fldName.c_str(), vtkDataObject::CELL)
                 || input->IsFieldNeeded(fldName.c_str(), vtkDataObject::POINT)
                )
                {
                    requestedFields.insert(fldName);
                }
            }

            dataq.set(channel, requestedFields);
        }
        else
        {
            dataq.clear(channel);
        }
    }

    return dataq.size();
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const catalyst::timeQuery& when)
{
    os << "Time = " << when.timeValue << ", index: " << when.timeIndex;
    return os;
}


// ************************************************************************* //

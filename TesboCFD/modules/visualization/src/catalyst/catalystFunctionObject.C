/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
    Copyright (C) 2018 CINECA
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

#include "catalystFunctionObject.H"
#include "catalystCoprocess.H"
#include "stringOps.H"
#include "OSHA1stream.H"
#include "OSspecific.H"  // For isFile, mkDir
#include "sigFpe.H"
#include "addToRunTimeSelectionTable.H"

#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>
#include <vtkSmartPointer.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(catalystFunctionObject, 0);

    addRemovableToRunTimeSelectionTable
    (
        functionObject,
        catalystFunctionObject,
        dictionary
    );

} // End namespace functionObjects
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::functionObjects::catalystFunctionObject::expand
(
    List<string>& scripts,
    const dictionary& dict
)
{
    label nscript = 0;

    forAll(scripts, scripti)
    {
        string& s = scripts[scripti];

        stringOps::inplaceExpand(s, dict, true, true);
        fileName::clean(s);  // Remove trailing, repeated slashes etc.

        if (Foam::isFile(s))
        {
            if (nscript != scripti)
            {
                scripts[nscript] = std::move(s);
            }
            ++nscript;
        }
    }

    const label nbad = (scripts.size() - nscript);

    if (nbad)
    {
        Warning
            << nl
            << "Could not resolve " << nbad << " of " << scripts.size()
            << " catalyst scripts - check your input" << nl << nl;
    }

    if (!nscript)
    {
        Warning
            << nl
            << "No catalyst scripts resolved - check your input" << nl << nl;
    }

    scripts.resize(nscript);

    return nscript;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::catalystFunctionObject::catalystFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    outputDir_("<case>/insitu"),
    scripts_(),
    adaptor_(),
    inputs_()
{
    if (postProcess)
    {
        // Disable for post-process mode.
        // Emit as FatalError for the try/catch in the caller.
        FatalError
            << type() << " disabled in post-process mode"
            << exit(FatalError);
    }

    if (!read(dict))
    {
        FatalError
            << type() << " with errors"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::catalystFunctionObject::~catalystFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::catalystFunctionObject::read(const dictionary& dict)
{
    functionObject::read(dict);

    const dictionary* inputDictPtr = dict.findDict("inputs");

    if (!inputDictPtr || inputDictPtr->empty())
    {
        inputs_.clear();

        WarningInFunction
            << type() << " missing inputs" << endl;

        return false;
    }

    int debugLevel = 0;
    if (dict.readIfPresent("debug", debugLevel))
    {
        catalyst::coprocess::debug = debugLevel;
    }

    if (UPstream::master())
    {
        fileName dir;
        if (dict.readIfPresent("mkdir", dir))
        {
            dir.expand();
            dir.clean();
            Foam::mkDir(dir);
        }
    }


    // Track changes in outputDir and/or scripts
    OSHA1stream osha1;
    osha1 << outputDir_ << scripts_;
    const SHA1Digest oldDigest(osha1.digest());

    dict.readIfPresent("outputDir", outputDir_);
    outputDir_.expand();
    outputDir_.clean();
    if (UPstream::master())
    {
        Foam::mkDir(outputDir_);
    }

    dict.readEntry("scripts", scripts_);    // Python scripts
    expand(scripts_, dict);                 // Expand and check availability


    // Any changes detected
    osha1.reset();
    osha1 << outputDir_ << scripts_;

    if (adaptor_ && oldDigest != osha1.digest())
    {
        adaptor_->reset(outputDir_, scripts_);
    }


    PtrList<catalyst::catalystInput> newList(inputDictPtr->size());

    label nInputs = 0;
    forAllConstIters(*inputDictPtr, iter)
    {
        const dictionary* subDictPtr = iter().dictPtr();
        if (iter().keyword().isPattern() || !subDictPtr)
        {
            continue;
        }

        auto input =
            catalyst::catalystInput::New
            (
                word(iter().keyword()),
                time_,
                *subDictPtr
            );

        // We may wish to perform additional validity or sanity checks on
        // the input before deciding to add it to the list.

        newList.set(nInputs, input);
        ++nInputs;
    }

    newList.resize(nInputs);
    inputs_.transfer(newList);

    if (inputs_.empty())
    {
        WarningInFunction
            << type() << " missing inputs" << endl;
    }
    else
    {
        Info<< type() << ' ' << name() << ':' << nl
            << "    output: " << outputDir_ << nl
            << "    scripts: " << scripts_ << nl
            << "    inputs:" << nl
            << "(" << nl;

        nInputs = 0;
        for (const auto& inp : inputs_)
        {
            if (nInputs++) Info<< nl;
            inp.print(Info);
        }

        Info<< ")" << nl;
    }

    return true;
}


bool Foam::functionObjects::catalystFunctionObject::execute()
{
    if (inputs_.empty())
    {
        return false;
    }

    // Disable any floating point trapping so that errors in
    // Catalyst pipelines do not kill the simulation
    // TODO: report that errors occurred?

    sigFpe::ignore sigFpeHandling; //<- disable in local scope

    if (!adaptor_)
    {
        adaptor_.reset(new catalyst::coprocess());
        adaptor_->reset(outputDir_, scripts_);
    }

    catalyst::dataQuery dataq(time_);

    label nChannels = 0;
    for (auto& inp : inputs_)
    {
        nChannels += inp.addChannels(dataq);
    }

    if (nChannels)
    {
        if (catalyst::coprocess::debug > 1)
        {
            Pout<< type() << ": query catalyst for " << nChannels
                << " channels" << nl;
        }

        nChannels = adaptor_().query(dataq);
    }

    if (catalyst::coprocess::debug > 1)
    {
        Pout<< type() << ": query catalyst for " << nChannels
            << " channels" << nl;
    }
    else if (catalyst::coprocess::debug)
    {
        Info<< type() << ": expecting data for " << nChannels
            << " channels" << nl;
    }

    if (!nChannels)
    {
        return true;
    }

    catalyst::outputChannels outputs;

    if (catalyst::coprocess::debug > 1)
    {
        Pout<< type() << ": converting input" << nl;
    }

    for (auto& inp : inputs_)
    {
        inp.convert(dataq, outputs);
    }

    if (catalyst::coprocess::debug > 1)
    {
        Pout<< type() << ": sending data for " << outputs.size()
            << " outputs" << nl;
    }

    if (outputs.size())
    {
        if (log)
        {
            Info<< type() << ": send data (";
            for (const word& channelName : outputs.sortedToc())
            {
                Info<< ' ' << channelName;
            }
            Info<< " )" << nl;
        }

        adaptor_().process(dataq, outputs);
    }


    if (catalyst::coprocess::debug > 1)
    {
        Pout<< type() << ": done step" << nl;
    }


    // Instead of relying on the destructor, manually restore the previous
    // SIGFPE state.
    // This is only to avoid compiler complaints about unused variables.

    sigFpeHandling.restore();

    return true;
}


bool Foam::functionObjects::catalystFunctionObject::write()
{
    return true;
}


bool Foam::functionObjects::catalystFunctionObject::end()
{
    // Only here for extra feedback
    if (log && adaptor_)
    {
        Info<< type() << ": Disconnecting ParaView Catalyst..." << nl;
    }

    adaptor_.clear();
    inputs_.clear();

    return true;
}


void Foam::functionObjects::catalystFunctionObject::updateMesh
(
    const mapPolyMesh&
)
{
    for (auto& inp : inputs_)
    {
        inp.update(polyMesh::TOPO_CHANGE);
    }
}


void Foam::functionObjects::catalystFunctionObject::movePoints
(
    const polyMesh&
)
{
    for (auto& inp : inputs_)
    {
        inp.update(polyMesh::POINTS_MOVED);
    }
}


// ************************************************************************* //

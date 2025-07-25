/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "adiosReader.H"
#include "adiosCore.H"

#include "dictionary.H"
#include "stringOps.H"
#include "IOstreams.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Output for map, unsorted map
    template<class K, class V>
    Ostream& operator<<(Ostream& os, const std::map<K, V>& tbl)
    {
        const label len = tbl.size();

        if (len)
        {
            // Size and start list delimiter
            os << nl << len << nl << token::BEGIN_LIST << nl;

            // Contents
            forAllConstIters(tbl, iter)
            {
                os << iter->first << ' ' << iter->second << nl;
            }

            os << token::END_LIST;    // End list delimiter
        }
        else
        {
            // Empty hash table
            os << len << token::BEGIN_LIST << token::END_LIST;
        }

        os.check(FUNCTION_NAME);
        return os;
    }


    // Output for Variable
    template<class T>
    Ostream& operator<<(Ostream& os, const adios2::Variable<T>& var)
    {
        if (var)
        {
            os.beginBlock();
            os.writeEntry("name",  var.Name());
            os.writeEntry("type",  var.Type());
            os.writeEntry("sizeof", var.Sizeof());

            os.writeEntry("steps", var.Steps());
            os.writeEntry("block", var.BlockID());

            os.writeEntry("selection", var.SelectionSize());

            os.writeEntry("nbytes", var.Sizeof()*var.SelectionSize());

            adios2::Dims count = var.Count();
            os << "count: " << count << nl;

            os.endBlock();
        }

        return os;
    }

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adiosFoam::adiosReader::reset()
{
    availableAttr_.clear();
    availableVars_.clear();

    labelByteSize_ = sizeof(label);
    scalarByteSize_ = sizeof(scalar);

    regionNames_.clear();
    clouds_.clear();
    fields_.clear();
}


void Foam::adiosFoam::adiosReader::scan(bool verbose)
{
    reset();

    if (!isGood())
    {
        if (verbose)
        {
            Info<< "adios-file not open" << endl;
        }
        return;
    }

    availableAttr_ = readIOPtr_->AvailableAttributes();
    availableVars_ = readIOPtr_->AvailableVariables();

    if (verbose)
    {
        Info<< "adios-file has "
            << availableVars_.size() << " variables and "
            << availableAttr_.size() << " attributes"
            << endl;
    }


    // Sizes of label, scalar
    label ival;
    if (readIntAttributeIfPresent(adiosFoam::foamAttribute/"label", ival))
    {
        labelByteSize_ = ival;
    }
    if (readIntAttributeIfPresent(adiosFoam::foamAttribute/"scalar", ival))
    {
        scalarByteSize_ = ival;
    }

    // Mandatory: /openfoam/regions
    regionNames_ = getStringListAttribute<word>
    (
        adiosFoam::foamAttribute/"regions"
    );


    for (const word& regName : regionNames_)
    {
        // Fields
        {
            const string prefixField(adiosFoam::fieldPath(regName));

            auto& tbl = fields_(regName);

            forAllConstIters(availableVars_, iter)
            {
                const fileName varName(iter->first);

                if (varName.ends_with("_0"))
                {
                    // The "_0" fields are addressed separately
                    continue;
                }

                const auto nameCmpts = stringOps::split(varName, '/');

                if (nameCmpts.size() == 3 && varName.startsWith(prefixField))
                {
                    // Match <regionName>/field/<fieldName>

                    tbl.insert
                    (
                        adiosFoam::fieldInfo
                        (
                            varName,
                            getStringAttribute(varName/"class")
                        )
                    );
                }
            }
        }


        // Clouds (optional): regionName/clouds
        {
            wordList cloudNames;
            readStringListAttributeIfPresent(regName/"clouds", cloudNames);

            auto& tbl = clouds_(regName);

            for (const word& cldName : cloudNames)
            {
                adiosFoam::cloudInfo item(regName, cldName, "");

                item.type() = getStringAttribute(item.fullName() / "class");

                tbl.insert(std::move(item));
            }
        }
    }


    // TODO: Parallel combine?

    if (verbose)
    {
        Info<< "Regions: " << flatOutput(regionNames_) << nl
            << "Fields: " << fields_ << nl
            << "Clouds: " << clouds_ << nl;
    }


    // mandatory: /openfoam/nRegions
    // - but could also just rely on region names
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adiosFoam::adiosReader::adiosReader(adios2::ADIOS& adios)
:
    adios_(adios),
    readIOPtr_(nullptr),
    readFilePtr_(nullptr),
    availableAttr_(),
    availableVars_(),
    labelByteSize_(sizeof(label)),
    scalarByteSize_(sizeof(scalar)),
    regionNames_(),
    clouds_(),
    fields_()
{}


Foam::adiosFoam::adiosReader::adiosReader(adios2::ADIOS& adios, const fileName& dataFile)
:
    adiosReader(adios)
{
    open(dataFile);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adiosFoam::adiosReader::~adiosReader()
{
    close();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::adiosFoam::adiosReader::open(const fileName& dataFile)
{
    // Close anything already open
    close();

    if (!readIOPtr_)
    {
        try
        {
            (void) adios_.DeclareIO("read");
        }
        catch (...)
        {
            // Need catch when in debug mode?
        }

        readIOPtr_.reset(new adios2::IO(adios_.AtIO("read")));
    }

    if (!readFilePtr_)
    {
        // Default Engine "BPFile"

        readFilePtr_.reset
        (
            new adios2::Engine
            (
                readIOPtr_->Open(dataFile, adios2::Mode::Read)
            )
            // inIO.SetEngine("BP3");
        );
    }

    // Error handling

    if (!readFilePtr_)
    {
        return false;
    }

    // Populate with general information
    scan(true);

    return isGood();
}


void Foam::adiosFoam::adiosReader::close()
{
    reset();

    if (readFilePtr_)
    {
        readFilePtr_->Close();
        readFilePtr_ = nullptr;
    }
}


bool Foam::adiosFoam::adiosReader::isGood() const
{
    return bool(readFilePtr_);
}


// ************************************************************************* //

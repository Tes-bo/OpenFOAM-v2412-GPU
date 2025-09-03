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

#include <fstream>
#include "foamVtmWriter.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Local Class * * * * * * * * * * * * * * * //

void Foam::vtk::vtmWriter::vtmEntry::clear()
{
    type_ = NONE;
    name_.clear();
    file_.clear();
}


bool Foam::vtk::vtmWriter::vtmEntry::good() const noexcept
{
    return
    (
        type_ == vtmEntry::BEGIN_BLOCK
     || type_ == vtmEntry::END_BLOCK
     || (type_ == vtmEntry::DATA && file_.size())
    );
}


bool Foam::vtk::vtmWriter::vtmEntry::write(vtk::formatter& format) const
{
    if (type_ == vtmEntry::BEGIN_BLOCK)
    {
        format.openTag(vtk::fileTag::BLOCK);
        if (name_.size())
        {
            format.xmlAttr("name", name_);
        }
        format.closeTag();

        return true;
    }
    else if (type_ == vtmEntry::END_BLOCK)
    {
        format.endBlock();
        return true;
    }
    else if (type_ == vtmEntry::DATA && file_.size())
    {
        format.openTag(vtk::fileTag::DATA_SET);

        if (name_.size())
        {
            format.xmlAttr("name", name_);
        }

        format.xmlAttr("file", file_);

        format.closeTag(true);  // Empty tag. ie, <DataSet ... />
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::vtk::vtmWriter::pruneEmpty()
{
    const label nEntries = entries_.size();

    label dst=0;

    for (label src=0; src < nEntries; ++src)
    {
        if (entries_[src].good())
        {
            if (dst != src)
            {
                entries_[dst] = std::move(entries_[src]);
            }
            ++dst;
        }
    }

    const bool changed = (dst != nEntries);
    entries_.resize(dst);

    return changed;
}


bool Foam::vtk::vtmWriter::pruneEmptyBlocks()
{
    bool pruned = false;

    const label nEntries = entries_.size();

    while (true)
    {
        bool changed = false;

        for (label i=0; i < nEntries; ++i)
        {
            vtmEntry& e = entries_[i];

            if (e.isType(vtmEntry::BEGIN_BLOCK))
            {
                for (label j=i+1; j < nEntries; ++j)
                {
                    if (entries_[j].isType(vtmEntry::END_BLOCK))
                    {
                        e.clear();
                        entries_[j].clear();

                        changed = true;
                        break;
                    }
                    else if (!entries_[j].isType(vtmEntry::NONE))
                    {
                        break;
                    }
                }
            }
        }

        if (changed)
        {
            pruned = true;
        }
        else
        {
            break;
        }
    }

    // Collapse single-entry blocks when the names allow it

    // Transcribe, removing NONE entries
    pruneEmpty();

    return pruned;
}


bool Foam::vtk::vtmWriter::collapseBlocks()
{
    bool collapsed = false;

    const label nEntries = entries_.size();

    for (label i=0; i < nEntries-2; ++i)
    {
        vtmEntry& b = entries_[i];    // begin
        vtmEntry& d = entries_[i+1];  // data
        vtmEntry& e = entries_[i+2];  // end

        if
        (
            b.isType(vtmEntry::BEGIN_BLOCK)
         && e.isType(vtmEntry::END_BLOCK)
         && d.isType(vtmEntry::DATA)
         && (d.name_.empty() || d.name_ == b.name_)
        )
        {
            d.name_ = std::move(b.name_);

            b.clear();
            e.clear();

            collapsed = true;
        }
    }

    pruneEmpty();

    return collapsed;
}


void Foam::vtk::vtmWriter::repair(bool collapse)
{
    // Add or remove END_BLOCK

    label depth = 0;
    label nEntries = 0;

    for (vtmEntry& e : entries_)
    {
        if (e.isType(vtmEntry::BEGIN_BLOCK))
        {
            ++depth;
        }
        else if (e.isType(vtmEntry::END_BLOCK))
        {
            --depth;

            if (depth < 0)
            {
                // Truncate now and exit
                entries_.resize(nEntries);
                break;
            }
        }
        else if (e.isType(vtmEntry::DATA))
        {
            if (e.file_.empty())
            {
                // Bad entry - reset to NONE
                e.clear();
            }
        }

        ++nEntries;
    }

    // Close any dangling blocks
    while (depth--)
    {
        entries_.append(vtmEntry::endblock());
    }

    blocks_.clear();
    pruneEmpty();

    if (collapse)
    {
        pruneEmptyBlocks();
        collapseBlocks();
    }
}


void Foam::vtk::vtmWriter::dump(Ostream& os) const
{
    label depth = 0;

    // Output format is a mix of dictionary and JSON
    // the only purpose being for diagnostics

    for (const vtmEntry& e : entries_)
    {
        switch (e.type_)
        {
            case vtmEntry::NONE:
            {
                os.indent();
                os  << "none" << nl;
                break;
            }
            case vtmEntry::DATA:
            {
                os.indent();
                os  << "{ \"name\" : " << e.name_
                    << ", \"file\" : " << e.file_ << " }" << nl;
                break;
            }
            case vtmEntry::BEGIN_BLOCK:
            {
                ++depth;
                os.beginBlock(e.name_);
                break;
            }
            case vtmEntry::END_BLOCK:
            {
                --depth;
                os.endBlock();
                os << nl;
                break;
            }
        }
    }

    for (label i=0; i < depth; ++i)
    {
        os.decrIndent();
    }

    if (depth > 0)
    {
        os << "# Had " << depth << " unclosed blocks" << nl;
    }
    if (depth < 0)
    {
        os << "# Had " << (-depth) << " too many end blocks" << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::vtmWriter::vtmWriter()
:
    vtmWriter(true)
{}


Foam::vtk::vtmWriter::vtmWriter(bool autoName)
:
    autoName_(autoName),
    hasTime_(false),
    entries_(),
    blocks_(),
    timeValue_(Zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::vtmWriter::clear()
{
    entries_.clear();
    blocks_.clear();

    timeValue_ = Zero;
    hasTime_ = false;
}


bool Foam::vtk::vtmWriter::empty() const
{
    for (const auto& e : entries_)
    {
        if (e.isType(vtmEntry::DATA) && e.name_.size())
        {
            return false;
        }
    }

    return true;
}


Foam::label Foam::vtk::vtmWriter::size() const
{
    label ndata = 0;

    for (const auto& e : entries_)
    {
        if (e.isType(vtmEntry::DATA) && e.file_.size())
        {
            ++ndata;
        }
    }

    return ndata;
}


void Foam::vtk::vtmWriter::setTime(scalar timeValue)
{
    timeValue_ = timeValue;
    hasTime_ = true;
}


void Foam::vtk::vtmWriter::setTime(const Time& t)
{
    timeValue_ = t.value();
    hasTime_ = true;
}


Foam::label Foam::vtk::vtmWriter::beginBlock(const word& blockName)
{
    entries_.append(vtmEntry::block(blockName));
    blocks_.append(blockName);

    return blocks_.size();
}


Foam::label Foam::vtk::vtmWriter::endBlock(const word& blockName)
{
    if (!blocks_.empty())
    {
        // Verify expected end tag
        if (!blockName.empty() && blockName != blocks_.back())
        {
            WarningInFunction
                << "expecting to end block '" << blockName
                << "' but found '" << blocks_.back() << "' instead"
                << endl;
        }

        blocks_.pop_back();
        entries_.push_back(vtmEntry::endblock());
    }

    return blocks_.size();
}


bool Foam::vtk::vtmWriter::append(const fileName& file)
{
    if (autoName_)
    {
        return append(fileName::stem(file), file);
    }

    return append(word::null, file);
}


bool Foam::vtk::vtmWriter::append
(
    const fileName& file,
    vtk::fileTag contentType
)
{
    if (autoName_)
    {
        return append(fileName::stem(file), file, contentType);
    }

    return append(word::null, file, contentType);
}


bool Foam::vtk::vtmWriter::append
(
    const word& name,
    const fileName& file
)
{
    if (file.empty())
    {
        return false;
    }

    entries_.append(vtmEntry::entry(name, file));
    return true;
}


bool Foam::vtk::vtmWriter::append
(
    const word& name,
    const fileName& file,
    vtk::fileTag contentType
)
{
    if (file.empty())
    {
        return false;
    }

    if (file.has_ext(vtk::fileExtension[contentType]))
    {
        entries_.append(vtmEntry::entry(name, file));
    }
    else
    {
        entries_.append
        (
            vtmEntry::entry
            (
                name,
                file + "." + vtk::fileExtension[contentType]
            )
        );
    }

    return true;
}


void Foam::vtk::vtmWriter::add
(
    const word& blockName,
    const fileName& prefix,
    const vtmWriter& other
)
{
    // Standard sanity repair (block ending), prune empty entries
    repair();

    beginBlock(blockName);

    label depth = 0;
    bool good = true;

    for (const vtmEntry& e : other.entries_)
    {
        switch (e.type_)
        {
            case vtmEntry::NONE:
            {
                break;
            }
            case vtmEntry::DATA:
            {
                if (e.good())
                {
                    entries_.append(e);

                    if (prefix.size())
                    {
                        fileName& f = entries_.last().file_;

                        f = prefix/f;
                    }
                }

                break;
            }
            case vtmEntry::BEGIN_BLOCK:
            {
                ++depth;
                entries_.append(e);
                break;
            }
            case vtmEntry::END_BLOCK:
            {
                good = (depth > 0);
                --depth;
                if (good)
                {
                    entries_.append(e);
                }
                break;
            }
        }

        if (!good) break;
    }

    while (depth--)
    {
        entries_.append(vtmEntry::endblock());
    }

    entries_.append(vtmEntry::endblock());

    if (!hasTime_ && other.hasTime_)
    {
        hasTime_ = true;
        timeValue_ = other.timeValue_;
    }
}


void Foam::vtk::vtmWriter::add
(
    const word& blockName,
    const vtmWriter& other
)
{
    add(blockName, fileName::null, other);
}


Foam::label Foam::vtk::vtmWriter::write(const fileName& file)
{
    std::ofstream os_;

    mkDir(file.path());

    if (file.has_ext(ext()))
    {
        os_.open(file);
    }
    else
    {
        os_.open(file + "." + ext());
    }

    auto format = vtk::newFormatter(os_, formatType::INLINE_ASCII);


    // Contents Header
    {
        format().xmlHeader();

        if (hasTime_)
        {
            format().xmlComment
            (
                "time='" + Foam::name(timeValue_) + "'"
            );
        }

        format().beginVTKFile<vtk::fileTag::MULTI_BLOCK>();
    }


    // Walk the block and dataset contents

    label depth = 0;
    label ndata = 0;

    for (const vtmEntry& e : entries_)
    {
        switch (e.type_)
        {
            case vtmEntry::DATA:
            {
                if (e.file_.empty())
                {
                    continue;  // Empty dataset is junk - skip
                }
                ++ndata;
                break;
            }
            case vtmEntry::BEGIN_BLOCK:
            {
                ++depth;
                break;
            }
            case vtmEntry::END_BLOCK:
            {
                --depth;
                break;
            }
            default:
            {
                continue;
                break;
            }
        }

        if (depth < 0)
        {
            // Too many end blocks - stop output now. Should we warn?
            break;
        }
        e.write(format());
    }

    // Close any dangling blocks
    while (depth--)
    {
        format().endBlock();
    }

    format().endTag(vtk::fileTag::MULTI_BLOCK);


    // FieldData for TimeValue
    if (hasTime_)
    {
        format()
            .beginFieldData()
            .writeTimeValue(timeValue_)
            .endFieldData();
    }

    format().endVTKFile();

    format.clear();
    os_.close();

    return ndata;
}


// ************************************************************************* //

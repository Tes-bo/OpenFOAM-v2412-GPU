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

#include "foamVtkSeriesWriter.H"
#include "Fstream.H"
#include "ListOps.H"
#include "stringOpsSort.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Get any single token.
    static inline bool getToken(ISstream& is, token& tok)
    {
        return (!is.read(tok).bad() && tok.good());
    }

    // Get two tokens.
    // The first one must be a ':' token, the second one is any value
    //
    // This corrsponds to the JSON  "key" : value syntax,
    // we trigger after reading the "key".
    static inline bool getValueToken(ISstream& is, token& tok)
    {
        return
        (
            // Token 1 = ':' separator
            (getToken(is, tok) && tok.isPunctuation(token::COLON))

            // Token 2 is the value
         && getToken(is, tok)
        );
    }


    // Sorting for fileNameInstant
    //   1. sort by value (time)
    //   2. natural sort (name)
    struct seriesLess
    {
        bool operator()(const fileNameInstant a, const fileNameInstant b) const
        {
            scalar val = compareOp<scalar>()(a.value(), b.value());
            if (val == 0)
            {
                return
                    stringOps::natural_sort::compare(a.name(), b.name()) < 0;
            }
            return val < 0;
        }
    };


    // Check if value is less than upper, with some tolerance.
    static inline bool lessThan(const scalar& val, const scalar& upper)
    {
        return (val < upper && Foam::mag(val - upper) > ROOTVSMALL);
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::fileName Foam::vtk::seriesWriter::base
(
    const fileName& outputName,
    char sep
)
{
    const auto dash = outputName.rfind(sep);

    // No separator? Or separator in path() instead of name()?
    if
    (
        std::string::npos == dash
     || std::string::npos != outputName.find('/', dash)
    )
    {
        // Warn?
        return outputName;
    }

    const auto dot = outputName.find('.', dash);

    if (std::string::npos == dot)
    {
        return outputName.substr(0, dash);
    }

    return outputName.substr(0, dash) + outputName.substr(dot);
}


Foam::word Foam::vtk::seriesWriter::suffix
(
    const fileName& file,
    char sep
)
{
    const auto dash = file.rfind(sep);

    // No separator? Or separator in path() instead of name()?
    if
    (
        std::string::npos == dash
     || std::string::npos != file.find('/', dash)
    )
    {
        // Warn?
        return "";
    }

    const auto dot = file.find('.', dash);

    if (std::string::npos == dot)
    {
        return file.substr(dash);
    }

    return file.substr(dash, (dot-dash));
}


Foam::Ostream& Foam::vtk::seriesWriter::print
(
    Ostream& os,
    const fileName& base,
    const UList<instant>& series,
    const char sep
)
{
    // Split the base into (stem, ext) components
    //
    // base = "path/file.vtm"
    //
    // stem = "file"
    // ext = ".vtm"

    const word stem = base.stem();
    const word ext = "." + base.ext();

    // Begin file-series (JSON)
    os  << "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n";

    // Track how many entries are remaining
    // - trailing commas on all but the final entry (JSON requirement)
    label nremain = series.size();

    // Each entry
    //   { "name" : "<stem><sep>name<ext>",  "time" : value }

    for (const instant& inst : series)
    {
        os  << "    { \"name\" : \""
            << stem << sep << inst.name() << ext
            << "\", \"time\" : " << inst.value() << " }";

        if (--nremain)
        {
            os  << ',';
        }
        os  << nl;
    }

    os  << "  ]\n}\n";

    return os;
}


Foam::Ostream& Foam::vtk::seriesWriter::print
(
    Ostream& os,
    const UList<fileNameInstant>& series
)
{
    // Begin file-series (JSON)
    os  << "{\n  \"file-series-version\" : \"1.0\",\n  \"files\" : [\n";

    // Track how many entries are remaining
    // - trailing commas on all but the final entry (JSON requirement)
    label nremain = series.size();

    // Each entry
    //   { "name" : "<file>",  "time" : <value> }

    for (const fileNameInstant& inst : series)
    {
        os  << "    { \"name\" : \""
            << inst.name().name()
            << "\", \"time\" : " << inst.value() << " }";

        if (--nremain)
        {
            os  << ',';
        }
        os  << nl;
    }

    os  << "  ]\n}\n";

    return os;
}


void Foam::vtk::seriesWriter::write
(
    const fileName& seriesName,
    const UList<instant>& series,
    const char sep
)
{
    mkDir(seriesName.path());

    autoPtr<OFstream> osPtr =
    (
        seriesName.has_ext("series")
      ? autoPtr<OFstream>::New(seriesName)
      : autoPtr<OFstream>::New(seriesName + ".series")
    );

    print(*osPtr, seriesName, series, sep);
}



void Foam::vtk::seriesWriter::write
(
    const fileName& seriesName,
    const UList<fileNameInstant>& series
)
{
    mkDir(seriesName.path());

    autoPtr<OFstream> osPtr =
    (
        seriesName.has_ext("series")
      ? autoPtr<OFstream>::New(seriesName)
      : autoPtr<OFstream>::New(seriesName + ".series")
    );

    print(*osPtr, series);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::vtk::seriesWriter::appendCheck(fileNameInstant inst)
{
    if (inst.name().empty())
    {
        return false;
    }

    const auto iter = existing_.find(inst.name());

    if (iter.good())
    {
        for (fileNameInstant& dst : entries_)
        {
            if (dst.name() == inst.name())
            {
                // Replace value
                dst.value() = inst.value();
                return true;
            }
        }
    }

    entries_.append(inst);
    existing_.insert(inst.name());

    return true;
}


bool Foam::vtk::seriesWriter::removeDuplicates()
{
    const label nElem = entries_.size();

    HashTable<label, fileName> filesSeen(2*nElem);

    bool changed = false;

    for (label elemi=0; elemi < nElem; ++elemi)
    {
        fileNameInstant& inst = entries_[elemi];

        if (inst.name().empty())
        {
            changed = true;
        }
        else
        {
            auto iter = filesSeen.find(inst.name());

            if (iter.good())
            {
                // Mark previous location as being superseded
                entries_[*iter].name().clear();
                changed = true;

                *iter = elemi;  // The latest with this name
            }
            else
            {
                filesSeen.insert(inst.name(), elemi);
            }
        }
    }


    if (changed)
    {
        label dsti = 0;
        for (label elemi=0; elemi < nElem; ++elemi)
        {
            fileNameInstant& src = entries_[elemi];

            if (!src.name().empty())
            {
                if (dsti != elemi)
                {
                    entries_[dsti] = std::move(src);
                }
                ++dsti;
            }
        }

        entries_.resize(dsti);
    }

    return (nElem != entries_.size());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::vtk::seriesWriter::load
(
    const fileName& seriesName,
    const bool checkFiles,
    const scalar restartTime
)
{
    clear();

    fileName seriesFile(seriesName);
    if (!seriesFile.has_ext("series"))
    {
        seriesFile.ext("series");
    }

    if (!isFile(seriesFile))
    {
        return size();
    }

    HashSet<fileName> filesOnDisk;

    if (checkFiles)
    {
        filesOnDisk.insert(Foam::readDir(seriesFile.path()));
    }


    // Parse JSON content:
    //
    // {
    //   "file-series-version" : "1.0",
    //   "files" : [
    //      { "name" : "abc", "time" : 123 },
    //      { "name" : "def", "time" : 345 }
    //    ]
    // }

    // Parsing states
    enum parse
    {
        NONE,         // Looking for "files"
        FILES_ARRAY,  // Saw "file" : '['
        ENTRY,        // Parsing in { "name" : ... }
        DONE,         // Saw a ']' while in FILES_ARRAY
        FAIL          // Something bad happened
    };

    // Track if "file" and "time" keys have been located
    unsigned instStatus = 0;
    fileNameInstant inst;

    token tok;

    IFstream is(seriesFile);

    for
    (
        parse state = parse::NONE;
        (state != parse::DONE && state != parse::FAIL)
     && getToken(is, tok);
        /*nil*/
    )
    {
        switch (state)
        {
            // Still scanning for initial "files" entry
            case parse::NONE :
            {
                if (tok.isString() && tok.stringToken() == "files")
                {
                    // Expect "files" : [ ...

                    if
                    (
                        getValueToken(is, tok)
                     && tok.isPunctuation(token::BEGIN_SQR)
                    )
                    {
                        state = parse::FILES_ARRAY;
                    }
                    else
                    {
                        state = parse::FAIL;
                    }
                }
            }
            break;

            // Parsing entries within "files" array
            case parse::FILES_ARRAY :
            {
                if (tok.isPunctuation())
                {
                    switch (tok.pToken())
                    {
                        // ',' - keep going (another entry)
                        case token::COMMA :
                            break;

                        // '{' - begin entry
                        case token::BEGIN_BLOCK :
                            state = parse::ENTRY;
                            instStatus = 0;
                            break;

                        // ']' - done array
                        case token::END_SQR :
                            state = parse::DONE;
                            break;

                        default:
                            state = parse::FAIL;
                            break;
                    }
                }
                else
                {
                    state = parse::FAIL;
                }
            }
            break;

            // Parsing an individual entry within "files"
            case parse::ENTRY :
            {
                if (tok.isPunctuation())
                {
                    switch (tok.pToken())
                    {
                        // ',' - keep going (another key/value pair)
                        case token::COMMA :
                            break;

                        // '}'
                        case token::END_BLOCK :
                        {
                            // Verify instant was properly parsed and
                            // is also valid
                            if
                            (
                                instStatus == 0x03
                             && lessThan(inst.value(), restartTime)
                             &&
                                (
                                    checkFiles
                                  ? filesOnDisk.found(inst.name())
                                  : true
                                )
                            )
                            {
                                appendCheck(inst);
                            }

                            state = parse::FILES_ARRAY;
                            instStatus = 0;
                        }
                        break;

                        default:
                            state = parse::FAIL;
                            break;
                    }
                }
                else if (tok.isString())
                {
                    // Expect "key" : value

                    const string key(tok.stringToken());

                    if (getValueToken(is, tok))
                    {
                        if ("name" == key)
                        {
                            if (tok.isString())
                            {
                                inst.name() = tok.stringToken();
                                instStatus |= 0x01;
                            }
                            else
                            {
                                state = parse::FAIL;
                            }
                        }
                        else if ("time" == key)
                        {
                            if (tok.isNumber())
                            {
                                inst.value() = tok.number();
                                instStatus |= 0x02;
                            }
                            else
                            {
                                state = parse::FAIL;
                            }
                        }
                    }
                    else
                    {
                        state = parse::FAIL;
                    }
                }
                else
                {
                    state = parse::FAIL;
                }
            }
            break;

            default:
                break;
        }
    }

    return size();
}


Foam::label Foam::vtk::seriesWriter::scan
(
    const fileName& seriesName,
    const scalar restartTime
)
{
    clear();

    const fileName path = seriesName.path();

    if (!isDir(path))
    {
        return size();
    }

    fileName seriesFile(seriesName);

    if (seriesName.has_ext("series"))
    {
        seriesFile.remove_ext();
    }

    const word stem = seriesFile.stem();
    const word ext = seriesFile.ext();

    // Accept "fileN.ext", "fileNN.ext", but reject "file.ext"
    const auto minLen = stem.length() + ext.length() + 1;

    const auto acceptName =
        [=](const fileName& file) -> bool
        {
            return
            (
                minLen < file.length()
             && file.has_ext(ext) && file.starts_with(stem)
            );
        };


    fileNameList files = subsetList(Foam::readDir(path), acceptName);

    // Names sorted so warnings appear less random
    Foam::sort(files, stringOps::natural_sort());

    // Scratch space for reading some of the file
    std::string header;

    scalar timeValue;

    bool warnings = false;

    for (const fileName& file : files)
    {
        std::ifstream is(path/file);

        if (!is)
        {
            continue;
        }

        // Read directly into the string
        // 1024 (12 lines of 80 chars) is plenty for all comments

        header.resize(1024);
        is.read(&(header.front()), header.size());
        const std::streamsize gcount = is.gcount();
        header.erase(gcount <= 0 ? 0 : gcount);

        // DebugInfo
        //     << "got header:\n=====\n" << header << "\n=====\n" << nl;


        // Look for time="...", time='...', or even time=... attribute

        auto begAttr = header.find("time=");

        if (string::npos == begAttr)
        {
            if (!warnings)
            {
                Info<< "No 'time=' comment attribute found:\n(" << nl;
                warnings = true;
            }
            Info<< "    " << file << nl;
            continue;
        }

        // Skip past the 'time='
        begAttr += 5;
        const char quote = header[begAttr];

        // Info<< "have time=" << int(begAttr) << nl;

        auto endAttr =
        (
            (quote == '"' || quote == '\'')
          ?
            // Quoted
            header.find(quote, ++begAttr)
          :
            // Unquoted
            header.find_first_of("\t\n\v\f\r ", begAttr)
        );


        if
        (
            string::npos != endAttr && begAttr < endAttr
         && readScalar
            (
                header.substr(begAttr, endAttr-begAttr),
                timeValue
            )
         && lessThan(timeValue, restartTime)
        )
        {
            // Success
            append(timeValue, file);
        }
    }

    if (warnings)
    {
        Info<< ")" << nl << nl;
    }

    // Don't trust the order. Sort by time and name instead.
    this->sort();

    return size();
}


bool Foam::vtk::seriesWriter::removeNewer(const scalar timeValue)
{
    // Rebuild hash as side-effect
    existing_.clear();

    label dsti = 0;

    const label nElem = entries_.size();

    for (label elemi=0; elemi < nElem; ++elemi)
    {
        fileNameInstant& src = entries_[elemi];

        if (!src.name().empty() && lessThan(src.value(), timeValue))
        {
            if (dsti != elemi)
            {
                entries_[dsti] = std::move(src);
                existing_.insert(entries_[dsti].name());
            }
            ++dsti;
        }
    }

    entries_.resize(dsti);

    return (nElem != entries_.size());
}


void Foam::vtk::seriesWriter::sort()
{
    Foam::sort(entries_, seriesLess());
}


// ************************************************************************* //

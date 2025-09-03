/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "ensightSurfaceReader.H"
#include "ensightCase.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ensightSurfaceReader, 0);
    addToRunTimeSelectionTable(surfaceReader, ensightSurfaceReader, fileName);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Extract timeset and fileset from split line information
// when the minElements has been satisfied.
// For example,
// ----
//     model: 1   some-geometry
//     model: 1 1 some-geometry
// ----
// would be split *after* the 'model:' resulting in these sub-strings:
//
//    ("1" "some-geometry")
//    ("1" "1" some-geometry")
//
// thus call extractTimeset with minElements == 2
//
template<class StringType>
static inline labelPair extractTimeset
(
    const SubStrings<StringType>& split,
    const std::size_t minElements
)
{
    ISpanStream is;

    labelPair result(-1, -1);
    if (split.size() >= minElements)
    {
        is.reset(split[0]);
        is >> result.first();

        if (split.size() > minElements)
        {
            is.reset(split[1]);
            is >> result.second();
        }
    }

    return result;
}

} // End namespace Foam


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

bool Foam::ensightSurfaceReader::readLine(ISstream& is, std::string& line)
{
    do
    {
        is.getLine(line);

        // Trim out any '#' comments (trailing or otherwise)
        const auto pos = line.find('#');
        if (pos != std::string::npos)
        {
            line.erase(pos);
        }
        stringOps::inplaceTrimRight(line);
    }
    while (line.empty() && is.good());

    return !line.empty();
}


void Foam::ensightSurfaceReader::checkSection
(
    const word& expected,
    const string& buffer,
    const ISstream& is
)
{
    // Be more generous with our expectations.
    // Eg, ensight specifies the "VARIABLE" entry,
    // but accepts "VARIABLES" as well.

    if (!expected.empty() && !buffer.starts_with(expected))
    {
        FatalIOErrorInFunction(is)
            << "Expected section header '" << expected
            << "' but read " << buffer << nl
            << exit(FatalIOError);
    }

    DebugInfo
        << "Read section header: " << buffer.c_str() << nl;
}


void Foam::ensightSurfaceReader::debugSection
(
    const word& expected,
    ISstream& is
)
{
    string buffer;
    readLine(is, buffer);

    checkSection(expected, buffer, is);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::ensightSurfaceReader::idTypes>
Foam::ensightSurfaceReader::readGeometryHeader
(
    ensightReadFile& is
) const
{
    string buffer;

    Pair<idTypes> idHandling(idTypes::NONE, idTypes::NONE);

    // Ensight Geometry File
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    // Description - 1
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    // "node id (off|assign|given|ignore)" - "given" is not actually supported
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    if (buffer.contains("ignore"))
    {
        idHandling.first() = idTypes::IGNORE;
    }
    else if (buffer.contains("given"))
    {
        idHandling.first() = idTypes::GIVEN;
    }

    // "element id (off|assign|given|ignore)"
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    if (buffer.contains("ignore"))
    {
        idHandling.second() = idTypes::IGNORE;
    }
    else if (buffer.contains("given"))
    {
        idHandling.second() = idTypes::GIVEN;
    }


    // "part" - but could also be an optional "extents"
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    if (buffer.contains("extents"))
    {
        // Optional extents - read and discard 6 floats
        // (xmin, xmax, ymin, ymax, zmin, zmax)

        is.skip<scalar>(6);

        // "part"
        is.read(buffer);
        DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;
    }

    // The part number
    label intValue;
    is.read(intValue);
    DebugInfo<< "part number: " << intValue << nl;

    // Part description / name
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    // "coordinates"
    is.read(buffer);
    DebugInfo<< "buffer [" << buffer.length() << "] " << buffer << nl;

    return idHandling;
}


void Foam::ensightSurfaceReader::readCase(ISstream& is)
{
    DebugInFunction << endl;

    enum ParseSection { UNKNOWN, FORMAT, GEOMETRY, VARIABLE, TIME, FILE };

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    string buffer;
    SubStrings<string> split;

    ParseSection parseState = ParseSection::UNKNOWN;

    // FORMAT
    // ~~~~~~

    debugSection("FORMAT", is);
    readLine(is, buffer);  // type: ensight gold
    parseState = ParseSection::FORMAT;

    // GEOMETRY
    // ~~~~~~~~
    debugSection("GEOMETRY", is);
    parseState = ParseSection::GEOMETRY;

    do
    {
        readLine(is, buffer);

        if (buffer.starts_with("VARIABLE"))
        {
            parseState = ParseSection::VARIABLE;
            break;
        }

        if (buffer.contains("change_coords_only"))
        {
            FatalIOErrorInFunction(is)
                << "No support for moving points, only topology change" << nl
                << exit(FatalIOError);
        }

        // Extract filename from GEOMETRY section
        // ====
        // model: [ts] [fs] filename  [change_coords_only [cstep]]
        //
        // ====

        // TBD:
        // check for "model:" vs "measured:" ?

        const auto pos_colon = buffer.find(':');

        if
        (
            (pos_colon == std::string::npos)
        )
        {
            FatalIOErrorInFunction(is)
                << "Error reading geometry 'model:'" << nl
                << exit(FatalIOError);
        }

        split = stringOps::splitSpace(buffer, pos_colon+1);

        if (split.empty())
        {
            FatalIOErrorInFunction(is)
                << "Error reading geometry 'model:'" << nl
                << exit(FatalIOError);
        }

        // With timeset? - need at least 2 entries
        meshTimeset_ = extractTimeset(split, 2);
        meshFileName_ = split[split.size()-1].str();

        DebugInfo << "mesh file:" << meshFileName_ << endl;
    }
    while (false);


    if (parseState != ParseSection::VARIABLE)
    {
        debugSection("VARIABLE", is);
        parseState = ParseSection::VARIABLE;
    }

    // Read the field description
    DynamicList<labelPair> dynFieldTimesets(16);
    DynamicList<word> dynFieldNames(16);
    DynamicList<string> dynFieldFileNames(16);

    // VARIABLE
    // ~~~~~~~~

    while (is.good())
    {
        readLine(is, buffer);

        if (buffer.starts_with("TIME"))
        {
            parseState = ParseSection::TIME;
            break;
        }

        // Read the field name and associated file name.
        // Eg,
        //     scalar per element: [ts] [fs]  p  data/********/p
        // but ignore
        //     scalar per node: [ts] [fs]  other  data/********/other

        const auto pos_colon = buffer.find(':');

        if (pos_colon == std::string::npos)
        {
            DebugInfo<< "ignore variable line: " << buffer << nl;
            continue;
        }

        // TODO? handle variable descriptions with spaces (they are quoted)

        split = stringOps::splitSpace(buffer, pos_colon+1);

        if (split.size() < 2)
        {
            WarningInFunction
                << "Error reading field file name, variable line: "
                << buffer << endl;
            continue;
        }

        auto pos_key = buffer.find("element");
        if ((pos_key == std::string::npos) || (pos_colon < pos_key))
        {
            DebugInfo<< "ignore variable line: " << buffer << nl;
            continue;
        }

        DebugInfo<< "variable line: " << buffer << nl;

        // With timeset? - need at least 3 entries
        dynFieldTimesets.push_back(extractTimeset(split, 3));

        dynFieldNames.push_back(split[split.size()-2].str());
        dynFieldFileNames.push_back(split[split.size()-1].str());
    }
    fieldTimesets_.transfer(dynFieldTimesets);
    fieldNames_.transfer(dynFieldNames);
    fieldFileNames_.transfer(dynFieldFileNames);

    DebugInfo
        << "fieldNames: " << fieldNames_ << nl
        << "fieldFileNames: " << fieldFileNames_ << nl;


    if (parseState != ParseSection::TIME)
    {
        FatalIOErrorInFunction(is)
            << "Did not find section header 'TIME'" << nl
            << exit(FatalIOError);
    }


    // Determine which unique timeset or fileset to expect

    labelHashSet expectTimeset;
    labelHashSet expectFileset;

    expectTimeset.insert(meshTimeset_.first());
    expectFileset.insert(meshTimeset_.second());

    for (const auto& tup : fieldTimesets_)
    {
        expectTimeset.insert(tup.first());
        expectFileset.insert(tup.second());
    }

    // Remove placeholders
    expectTimeset.erase(-1);
    expectFileset.erase(-1);


    DebugInfo
        << "expect timesets: " << flatOutput(expectTimeset) << nl
        << "expect filesets: " << flatOutput(expectFileset) << nl;


    // TIME
    // ~~~~
    // style 1:
    // ====
    // time set:              <int>  [description]
    // number of steps:       <int>
    // filename start number: <int>
    // filename increment:    <int>
    // time values: time_1 .. time_N
    // ====
    //
    // style 2:
    // ====
    // time set:              <int>  [description]
    // number of steps:       <int>
    // filename numbers:      int_1 .. int_N
    // time values: time_1 .. time_N
    // ====
    //
    // style 3:
    // ====
    // time set:              <int>  [description]
    // number of steps:       <int>
    // filename numbers file: <filename>
    // time values file:      <filename>
    // ====


    // Currently only handling style 1, style 2
    // and only a single time set

    // time set = 1
    {
        // time set: <int>
        {
            readLine(is, buffer);
        }

        // number of steps: <int>
        label nTimes = 0;
        {
            readLine(is, buffer);
            split = stringOps::splitSpace(buffer);
            readFrom(split.back(), nTimes);
        }

        // filename start number: <int>
        // filename increment: <int>
        //
        // OR:
        // filename numbers: ...

        readLine(is, buffer);
        auto pos_colon = buffer.find(':');

        if (buffer.contains("numbers:"))
        {
            // Split out trailing values...
            split = stringOps::splitSpace(buffer, pos_colon+1);

            fileNumbers_.resize_nocopy(nTimes);

            label numRead = 0;
            while (numRead < nTimes)
            {
                for (const auto& chunk : split)
                {
                    std::string str(chunk.str());

                    if (!Foam::readLabel(str, fileNumbers_[numRead]))
                    {
                        FatalIOErrorInFunction(is)
                            << "Could not parse label: " << str << nl
                            << exit(FatalIOError);
                    }
                    ++numRead;

                    if (numRead == nTimes)
                    {
                        break;
                    }
                }

                // Get more input
                if (numRead < nTimes)
                {
                    readLine(is, buffer);
                    split = stringOps::splitSpace(buffer);
                }
            }

            timeStartIndex_ = 0;
            timeIncrement_ = 0;
            fileNumbers_.resize(numRead);
        }
        else
        {
            // filename start number: <int>
            split = stringOps::splitSpace(buffer);
            readFrom(split.back(), timeStartIndex_);

            // filename increment: <int>
            readLine(is, buffer);
            split = stringOps::splitSpace(buffer);
            readFrom(split.back(), timeIncrement_);

            fileNumbers_.clear();
        }

        DebugInfo
            << "nTimes: " << nTimes
            << " start-index: " << timeStartIndex_
            << " increment: " << timeIncrement_
            << " file numbers: " << flatOutput(fileNumbers_) << nl;


        // time values: time_1 .. time_N
        readLine(is, buffer);

        // Split out trailing values...
        {
            const auto pos_colon = buffer.find(':');
            const auto pos_key = buffer.find("values");

            if
            (
                (pos_colon == std::string::npos)
             || (pos_key == std::string::npos) || (pos_colon < pos_key)
            )
            {
                split.clear();
            }
            else
            {
                split = stringOps::splitSpace(buffer, pos_colon+1);
            }
        }

        timeValues_.resize_nocopy(nTimes);

        label numRead = 0;
        while (numRead < nTimes)
        {
            for (const auto& chunk : split)
            {
                auto& inst = timeValues_[numRead];

                // Retain character representation
                inst.name() = word(chunk.str());

                if (!Foam::readScalar(inst.name(), inst.value()))
                {
                    FatalIOErrorInFunction(is)
                        << "Could not parse scalar: " << inst.name() << nl
                        << exit(FatalIOError);
                }
                ++numRead;

                if (numRead == nTimes)
                {
                    break;
                }
            }

            // Get more input
            if (numRead < nTimes)
            {
                readLine(is, buffer);
                split = stringOps::splitSpace(buffer);
            }
        }

        timeValues_.resize(numRead);
    }

    // Not yet:

    // FILE
    // ~~~~
    // file set:        <int>
    // filename index:  <int>   - file index number in the file name
    // number of steps: <int>
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightSurfaceReader::ensightSurfaceReader
(
    const fileName& fName,
    const dictionary& options
)
:
    surfaceReader(fName, options),
    masterOnly_
    (
        UPstream::parRun()
     && options.getOrDefault("masterOnly", false)
    ),
    readFormat_(IOstreamOption::ASCII),  // Placeholder value
    baseDir_(fName.path()),
    meshTimeset_(-1,-1),
    timeStartIndex_(0),
    timeIncrement_(1)
{
    if (options.getOrDefault("debug", false))
    {
        debug |= 1;
    }

    if (!masterOnly_ || UPstream::master(UPstream::worldComm))
    {
        IFstream is(fName);
        readCase(is);
    }

    if (masterOnly_ && UPstream::parRun())
    {
        Pstream::broadcasts
        (
            UPstream::worldComm,
            meshTimeset_,
            meshFileName_,
            fieldTimesets_,
            fieldNames_,
            fieldFileNames_,
            timeStartIndex_,
            timeIncrement_,
            fileNumbers_,
            timeValues_
        );
    }
}


// * * * * * * * * * * * * * Public Member Functions   * * * * * * * * * * * //

Foam::meshedSurface Foam::ensightSurfaceReader::readGeometry
(
    const fileName& geometryFile,
    const label timeIndex
)
{
    DebugInFunction << endl;

    {
        // Auto-detect ascii/binary format
        ensightReadFile is(geometryFile);

        // Format detected from the geometry
        readFormat_ = is.format();

        // For transient single-file
        is.seekTime(timeIndex);

        DebugInfo
            << "File: " << is.name()
            << " format: "
            << IOstreamOption::formatNames[readFormat_] << endl;

        Pair<idTypes> idHandling = readGeometryHeader(is);

        label nPoints;
        is.read(nPoints);

        DebugInfo
            << "nPoints: " << nPoints << nl;

        if (idHandling.first() == idTypes::GIVEN)
        {
            WarningInFunction
                << "Treating node id 'given' as being 'ignore'" << nl
                << "If something fails, this could be the reason" << nl
                << endl;

            idHandling.first() = idTypes::IGNORE;
        }

        if (idHandling.first() == idTypes::IGNORE)
        {
            DebugInfo
                << "Ignore " << nPoints << " node ids" << nl;

            // Read and discard labels
            is.skip<label>(nPoints);
        }

        pointField points;
        is.readPoints(nPoints, points);


        // Read faces - may be a mix of tria3, quad4, nsided
        DynamicList<face> dynFaces(nPoints/3);
        DynamicList<faceInfoTuple> faceTypeInfo(16);

        string buffer;

        while (is.good())
        {
            // The element type
            is.read(buffer);

            if (!is.good())
            {
                break;
            }
            else if (buffer.contains("BEGIN TIME STEP"))
            {
                // Graciously handle a miscued start
                continue;
            }
            else if (buffer.contains("END TIME STEP"))
            {
                // END TIME STEP is a valid means to terminate input
                break;
            }

            if
            (
                buffer
             == ensightFaces::elemNames[ensightFaces::elemType::TRIA3]
            )
            {
                label elemCount;
                is.read(elemCount);

                faceTypeInfo.emplace_back
                (
                    ensightFaces::elemType::TRIA3,
                    elemCount
                );

                DebugInfo
                    << "faceType <"
                    << ensightFaces::elemNames[ensightFaces::elemType::TRIA3]
                    << "> count: "
                    << elemCount << nl;

                if
                (
                    idHandling.second() == idTypes::IGNORE
                 || idHandling.second() == idTypes::GIVEN
                )
                {
                    DebugInfo
                        << "Ignore " << elemCount << " element ids" << nl;

                    // Read and discard labels
                    is.skip<label>(elemCount);
                }

                // Extend and fill the new trailing portion
                const label startElemi = dynFaces.size();
                dynFaces.resize(startElemi+elemCount, face(3));  // tria3
                faceList::subList myElements = dynFaces.slice(startElemi);

                for (auto& f : myElements)
                {
                    for (label& fp : f)
                    {
                        is.read(fp);
                    }
                }
            }
            else if
            (
                buffer
             == ensightFaces::elemNames[ensightFaces::elemType::QUAD4]
            )
            {
                label elemCount;
                is.read(elemCount);

                faceTypeInfo.emplace_back
                (
                    ensightFaces::elemType::QUAD4,
                    elemCount
                );

                DebugInfo
                    << "faceType <"
                    << ensightFaces::elemNames[ensightFaces::elemType::QUAD4]
                    << "> count: "
                    << elemCount << nl;

                if
                (
                    idHandling.second() == idTypes::IGNORE
                 || idHandling.second() == idTypes::GIVEN
                )
                {
                    DebugInfo
                        << "Ignore " << elemCount << " element ids" << nl;

                    // Read and discard labels
                    is.skip<label>(elemCount);
                }

                // Extend and fill the new trailing portion
                const label startElemi = dynFaces.size();
                dynFaces.resize(startElemi + elemCount, face(4));  // quad4
                faceList::subList myElements = dynFaces.slice(startElemi);

                for (auto& f : myElements)
                {
                    for (label& fp : f)
                    {
                        is.read(fp);
                    }
                }
            }
            else if
            (
                buffer
             == ensightFaces::elemNames[ensightFaces::elemType::NSIDED]
            )
            {
                label elemCount;
                is.read(elemCount);

                faceTypeInfo.emplace_back
                (
                    ensightFaces::elemType::NSIDED,
                    elemCount
                );

                DebugInfo
                    << "faceType <"
                    << ensightFaces::elemNames[ensightFaces::elemType::NSIDED]
                    << "> count: " << elemCount << nl;

                if
                (
                    idHandling.second() == idTypes::IGNORE
                 || idHandling.second() == idTypes::GIVEN
                )
                {
                    DebugInfo
                        << "Ignore " << elemCount << " element ids" << nl;

                    // Read and discard labels
                    is.skip<label>(elemCount);
                }

                // Extend and fill the new trailing portion
                const label startElemi = dynFaces.size();
                dynFaces.resize(startElemi + elemCount);
                faceList::subList myElements = dynFaces.slice(startElemi);

                for (auto& f : myElements)
                {
                    label nVerts;
                    is.read(nVerts);

                    f.resize(nVerts);
                }

                for (auto& f : myElements)
                {
                    for (label& fp : f)
                    {
                        is.read(fp);
                    }
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unknown face type: <" << buffer.c_str()
                        << ">. Stopping read and continuing with current "
                        << "elements only" << endl;
                }
                break;
            }
        }

        // From 1-based Ensight addressing to 0-based OF addressing
        for (auto& f : dynFaces)
        {
            for (label& fp : f)
            {
                --fp;
            }
        }

        faceTypeInfo_.transfer(faceTypeInfo);
        faceList faces(std::move(dynFaces));

        DebugInfo
            << "read nFaces: " << faces.size() << nl
            << "file schema: " << faceTypeInfo_ << nl;

        return meshedSurface(std::move(points), std::move(faces));
    }
}


const Foam::meshedSurface& Foam::ensightSurfaceReader::geometry
(
    const label timeIndex
)
{
    DebugInFunction << endl;

    if (!surfPtr_)
    {
        surfPtr_.reset(new meshedSurface);
        auto& surf = *surfPtr_;

        fileName geomFile
        (
            baseDir_
          / ensightCase::expand_mask(meshFileName_, timeIndex)
        );

        if (!masterOnly_ || UPstream::master(UPstream::worldComm))
        {
            surf = readGeometry(geomFile, timeIndex);
        }

        if (masterOnly_ && UPstream::parRun())
        {
            // Note: don't need faceTypeInfo_ on (non-reading) ranks
            Pstream::broadcast(surf, UPstream::worldComm);
        }
    }

    return *surfPtr_;
}


Foam::instantList Foam::ensightSurfaceReader::times() const
{
    return timeValues_;
}


Foam::wordList Foam::ensightSurfaceReader::fieldNames
(
    const label timeIndex
) const
{
    return fieldNames_;
}


Foam::tmp<Foam::Field<Foam::scalar>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const scalar& refValue
) const
{
    return readField<scalar>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::vector>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const vector& refValue
) const
{
    return readField<vector>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::sphericalTensor>>
Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const sphericalTensor& refValue
) const
{
    return readField<sphericalTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::symmTensor>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const symmTensor& refValue
) const
{
    return readField<symmTensor>(timeIndex, fieldIndex);
}


Foam::tmp<Foam::Field<Foam::tensor>> Foam::ensightSurfaceReader::field
(
    const label timeIndex,
    const label fieldIndex,
    const tensor& refValue
) const
{
    return readField<tensor>(timeIndex, fieldIndex);
}


// ************************************************************************* //

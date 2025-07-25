/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "ensightCase.H"
#include "Time.H"
#include "cloud.H"
#include "IOmanip.H"
#include "OSstream.H"
#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::ensightCase::dataDirName  = "data";
const char* Foam::ensightCase::geometryName = "geometry";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::ensightCase::mask(const int nwidth)
{
    if (nwidth < 1)
    {
        return word();
    }

    return word(std::string(nwidth, '*'), false);  // stripping=false
}


Foam::word Foam::ensightCase::padded(const int nwidth, const label index)
{
    if (nwidth < 1)
    {
        return Foam::name(index);
    }

    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(nwidth) << index;

    return word(oss.str(), false);  // stripping=false
}


void Foam::ensightCase::setTimeFormat
(
    OSstream& os,
    IOstreamOption::floatFormat timeFmt,
    const int timePrec
)
{
    os.setf(std::ios_base::left);
    os.setf
    (
        std::ios_base::fmtflags(timeFmt),
        std::ios_base::floatfield
    );

    if (timePrec > 0)
    {
        os.precision(timePrec);
    }
}


void Foam::ensightCase::setTimeFormat
(
    OSstream& os,
    const ensightCase::options& opts
)
{
    os.setf(std::ios_base::left);
    os.setf
    (
        std::ios_base::fmtflags(opts.timeFormat()),
        std::ios_base::floatfield
    );

    os.precision(opts.timePrecision());
}


void Foam::ensightCase::printTimeset
(
    OSstream& os,
    const label ts,
    const scalar timeValue
)
{
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << 1 << nl;

    // Single value - starts at index 0
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl
        << "time values:" << nl;

    os  << "    " << timeValue
        << nl << nl;
}


void Foam::ensightCase::printTimeset
(
    OSstream& os,
    const label ts,
    const UList<scalar>& values
)
{
    label pos_(0);

    os
        << "time set:               " << ts << nl
        << "number of steps:        " << values.size() << nl;

    // Assume contiguous numbering - starts at index 0
    os  << "filename start number:  0" << nl
        << "filename increment:     1" << nl;


    os  << "time values:" << nl;
    pos_ = 0;
    for (const scalar val : values)
    {
        if (pos_ == 6)
        {
            os  << nl;
            pos_ = 0;
        }
        ++pos_;

        os  << ' ' << setf(ios_base::right) << setw(12) << val;
    }
    os  << nl << nl;
}


void Foam::ensightCase::printTimeset
(
    OSstream& os,
    const label ts,
    const UList<scalar>& values,
    const bitSet& indices
)
{
    label pos_(0);

    // Check if continuous numbering can be used
    if
    (
        values.empty()
     || (indices.size() == values.size() && indices.all())
    )
    {
        // Can simply emit as 0-based with increment
        printTimeset(os, ts, values);
        return;
    }


    // Generate time set
    os
        << "time set:               " << ts << nl
        << "number of steps:        " << indices.count() << nl;


    os  << "filename numbers:" << nl;
    pos_ = 0;
    for (const label idx : indices)
    {
        if (pos_ == 6)
        {
            os  << nl;
            pos_ = 0;
        }
        ++pos_;

        os  << ' ' << setf(ios_base::right) << setw(8) << idx;
    }
    os  << nl;


    os  << "time values:" << nl;
    pos_ = 0;
    for (const label idx : indices)
    {
        if (pos_ == 6)
        {
            os  << nl;
            pos_ = 0;
        }
        ++pos_;

        os  << ' ' << setf(ios_base::right) << setw(12) << values[idx];
    }
    os  << nl << nl;
}


// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

Foam::fileName Foam::ensightCase::dataDir() const
{
    return ensightDir_/dataDirName;
}


void Foam::ensightCase::initialize()
{
    if (UPstream::master())
    {
        // EnSight and EnSight/data directories must exist

        // We may wish to retain old data
        // eg, convert new results or a particular time interval
        // OR remove everything

        if (Foam::isDir(ensightDir_))
        {
            if (options_->overwrite())
            {
                Foam::rmDir(ensightDir_);
            }
            else
            {
                DetailInfo
                    << "Warning: re-using existing directory" << nl
                    << "    " << ensightDir_ << endl;
            }
        }

        // Create ensight and data directories
        Foam::mkDir(dataDir());

        // The case file is always ASCII
        os_.reset(new OFstream(ensightDir_/caseName_, IOstreamOption::ASCII));
        ensightCase::setTimeFormat(*os_, *options_);  // Format options

        writeHeader();
    }
}


Foam::label Foam::ensightCase::checkTimeset(const labelHashSet& lookup) const
{
    // assume the worst
    label ts = -1;

    // work on a copy
    labelHashSet tsTimes(lookup);
    tsTimes.erase(-1);

    if (tsTimes.empty())
    {
        // no times needed
        ts = 0;
    }
    else if (tsTimes.size() == timesUsed_.size())
    {
        forAllConstIters(timesUsed_, iter)
        {
            tsTimes.erase(iter.key());
        }

        // OR
        // tsTimes.unset(timesUsed_.toc());

        if (tsTimes.empty())
        {
            ts = 1; // can use timeset 1
        }
    }

    return ts;
}


void Foam::ensightCase::writeHeader() const
{
    if (os_)  // True on master only
    {
        this->rewind();
        *os_
            << "FORMAT" << nl
            << "type: ensight gold" << nl;
    }
}


Foam::scalar Foam::ensightCase::writeTimeset() const
{
    const label ts = 1;

    const labelList indices(timesUsed_.sortedToc());
    label count = indices.size();

    // correct for negative starting values
    scalar timeCorrection = timesUsed_[indices[0]];
    if (timeCorrection < 0)
    {
        timeCorrection = -timeCorrection;
        Info<< "Correcting time values. Adding " << timeCorrection << endl;
    }
    else
    {
        timeCorrection = 0;
    }


    *os_
        << "time set:               " << ts << nl
        << "number of steps:        " << count << nl;

    if (indices[0] == 0 && indices[count-1] == count-1)
    {
        // looks to be contiguous numbering
        *os_
            << "filename start number:  " << 0 << nl
            << "filename increment:     " << 1 << nl;
    }
    else
    {
        *os_
            << "filename numbers:" << nl;

        count = 0;
        for (const label idx : indices)
        {
            *os_ << ' ' << setw(12) << idx;

            if (++count % 6 == 0)
            {
                *os_ << nl;
            }
        }

        if (count)
        {
            *os_ << nl;
        }
    }


    *os_ << "time values:" << nl;

    count = 0;
    for (const label idx : indices)
    {
        *os_ << ' ' << setw(12) << timesUsed_[idx] + timeCorrection;

        if (++count % 6 == 0)
        {
            *os_ << nl;
        }
    }
    if (count)
    {
        *os_ << nl;
    }

    return timeCorrection;
}


void Foam::ensightCase::writeTimeset
(
    const label ts,
    const labelHashSet& lookup,
    const scalar timeCorrection
) const
{
    // Make a copy
    labelHashSet hashed(lookup);
    hashed.erase(-1);

    const labelList indices(hashed.sortedToc());
    label count = indices.size();

    *os_
        << "time set:               " << ts << nl
        << "number of steps:        " << count  << nl
        << "filename numbers:" << nl;

    count = 0;
    for (const label idx : indices)
    {
        *os_ << ' ' << setw(12) << idx;

        if (++count % 6 == 0)
        {
            *os_ << nl;
        }
    }

    if (count)
    {
        *os_ << nl;
    }

    *os_ << "time values:" << nl;

    count = 0;
    for (const label idx : indices)
    {
        *os_ << ' ' << setw(12) << timesUsed_[idx] + timeCorrection;

        if (++count % 6 == 0)
        {
            *os_ << nl;
        }
    }
    if (count)
    {
        *os_ << nl;
    }
}


void Foam::ensightCase::noteGeometry(const bool moving) const
{
    if (moving)
    {
        geomTimes_.insert(timeIndex_);
    }
    else
    {
        geomTimes_.insert(-1);
    }

    changed_ = true;
}


void Foam::ensightCase::noteCloud(const word& cloudName) const
{
    // Force into existence
    if (!cloudVars_.found(cloudName))
    {
        cloudVars_.emplace(cloudName);
    }
    cloudTimes_.insert(timeIndex_);

    changed_ = true;
}


void Foam::ensightCase::noteCloud
(
    const word& cloudName,
    const word& varName,
    const char* ensightType
) const
{
    if (cloudVars_.found(cloudName))
    {
        if (cloudVars_[cloudName].insert(varName, ensightType))
        {
            changed_ = true;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Tried to add a cloud variable for writing"
            << " - without having added a cloud"
            << abort(FatalError);
    }
}


void Foam::ensightCase::noteVariable
(
    const word& varName,
    const char* ensightType
) const
{
    if (variables_.insert(varName, ensightType))
    {
        changed_ = true;
    }
}


Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::createDataFile
(
    const word& name
) const
{
    if (UPstream::master())
    {
        // The data/ITER subdirectory must exist
        // Note that data/ITER is indeed a valid ensight::FileName

        const fileName outdir = dataDir()/padded(timeIndex_);
        Foam::mkDir(outdir);

        return autoPtr<ensightFile>::New(outdir, name, format());
    }

    return nullptr;
}


Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::createCloudFile
(
    const word& cloudName,
    const word& name
) const
{
    if (UPstream::master())
    {
        // Write
        // eg -> "data/********/lagrangian/<cloudName>/positions"
        // or -> "lagrangian/<cloudName>/********/positions"
        // TODO? check that cloudName is a valid ensight filename
        const fileName outdir =
        (
            separateCloud()
          ? (ensightDir_ / cloud::prefix / cloudName / padded(timeIndex_))
          : (dataDir() / padded(timeIndex_) / cloud::prefix / cloudName)
        );

        Foam::mkDir(outdir); // should be unnecessary after newCloud()

        return autoPtr<ensightFile>::New(outdir, name, format());
    }

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightCase::ensightCase
(
    const fileName& ensightDir,
    const word& caseName,
    const ensightCase::options& opts
)
:
    options_(new options(opts)),
    os_(nullptr),
    ensightDir_(ensightDir),
    caseName_(caseName + ".case"),
    changed_(false),
    timeIndex_(0),
    timeValue_(0)
{
    initialize();
}


Foam::ensightCase::ensightCase
(
    const fileName& ensightDir,
    const word& caseName,
    const IOstreamOption::streamFormat fmt
)
:
    options_(new options(fmt)),
    os_(nullptr),
    ensightDir_(ensightDir),
    caseName_(caseName + ".case"),
    changed_(false),
    timeIndex_(0),
    timeValue_(0)
{
    initialize();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightCase::nextTime(const scalar value)
{
    // use next available index
    setTime(value, timesUsed_.size());
}


void Foam::ensightCase::nextTime(const instant& t)
{
    nextTime(t.value());
}


void Foam::ensightCase::setTime(const scalar value, const label index)
{
    timeIndex_ = index;
    timeValue_ = value;

    if (UPstream::master())
    {
        // The data/ITER subdirectory must exist
        // Note that data/ITER is indeed a valid ensight::FileName

        const fileName outdir = dataDir()/padded(timeIndex_);
        Foam::mkDir(outdir);

        // place a timestamp in the directory for future reference
        OFstream timeStamp(outdir/"time");
        timeStamp
            << "#  index  time" << nl
            << outdir.name() << ' ' << timeValue_ << nl;
    }

    // Record of time index/value used
    timesUsed_.set(index, value);
}


void Foam::ensightCase::setTime(const instant& t, const label index)
{
    setTime(t.value(), index);
}


void Foam::ensightCase::write() const
{
    if (!os_) return; // master only

    // geometry timeset
    const bool staticGeom = (geomTimes_.size() == 1 && geomTimes_.found(-1));
    label tsGeom = staticGeom ? 0 : checkTimeset(geomTimes_);

    // geometry index, when mesh is not moving but stored under data/XXX/
    label meshIndex = -1;

    // cloud timeset
    label tsCloud = checkTimeset(cloudTimes_);

    // Increment time-sets to the correct indices
    if (tsGeom < 0)
    {
        tsGeom = 2; // Next available timeset

        // Saved under data/XXX/geometry, but not actually moving
        if (geomTimes_.size() == 1)
        {
            tsGeom = 0;
            meshIndex = *(geomTimes_.begin());
        }
    }
    if (tsCloud < 0)
    {
        tsCloud = 1 + std::max(label(1), tsGeom);  // Next available timeset
    }

    writeHeader();


    // data mask: eg "data/******"
    const fileName dataMask = (dataDirName/mask());

    //
    // GEOMETRY
    //
    if (!geomTimes_.empty() || !cloudTimes_.empty())
    {
        // start of variables
        *os_
            << nl
            << "GEOMETRY" << nl;
    }

    if (staticGeom)
    {
        // Static mesh: store under data/constant/geometry
        *os_
            << setw(16)  << "model:"
            << (dataDirName/word("constant")/geometryName).c_str()
            << nl;
    }
    else if (meshIndex >= 0)
    {
        // Not really moving, but stored under data/XXXX/geometry
        *os_
            << setw(16)  << "model:"
            << (dataDirName/padded(meshIndex)/geometryName).c_str()
            << nl;
    }
    else if (!geomTimes_.empty())
    {
        // Moving
        *os_
            << word::printf("model: %-9d", tsGeom) // width 16 (no quotes)
            << (dataMask/geometryName).c_str()
            << nl;
    }

    // Clouds and cloud variables
    const wordList cloudNames(cloudVars_.sortedToc());

    for (const word& cloudName : cloudNames)
    {
        const fileName masked =
        (
            separateCloud()
          ? (cloud::prefix / cloudName / mask())
          : (dataMask / cloud::prefix / cloudName)
        );

        *os_
            << word::printf("measured: %-6d", tsCloud) // width 16 (no quotes)
            << (masked/"positions").c_str()
            << nl;
    }


    //
    // VARIABLE
    //
    if (variables_.size() || cloudVars_.size())
    {
        // Start of variables
        *os_
            << nl
            << "VARIABLE" << nl;
    }


    // Field variables (always use timeset 1)
    // NB: The output file name is stricter than the variable name

    for (const word& varName : variables_.sortedToc())
    {
        const string& ensType = variables_[varName];

        *os_
            << ensType.c_str()
            <<
            (
                (nodeVariables_.found(varName) || nodeValues())
              ? " per node:    1  "  // time-set 1
              : " per element: 1  "  // time-set 1
            )
            << setw(15) << varName << ' '
            << (dataMask/ensight::FileName(varName)).c_str() << nl;
    }


    // Clouds and cloud variables (using cloud timeset)
    // Write
    // as -> "data/********/lagrangian/<cloudName>/positions"
    // or -> "lagrangian/<cloudName>/********/positions"
    // NB: The output file name is stricter than the variable name

    label cloudNo = 0;
    for (const word& cloudName : cloudNames)
    {
        const fileName masked =
        (
            separateCloud()
          ? (cloud::prefix / cloudName / mask())
          : (dataMask / cloud::prefix / cloudName)
        );

        const HashTable<string>& vars = cloudVars_[cloudName];

        for (const word& varName : vars.sortedToc())
        {
            const string& ensType = vars[varName];

            // prefix variables with 'c' (cloud) and cloud index
            *os_
                << ensType.c_str() << " per "
                << word::printf("measured node: %-5d", tsCloud) // width 20
                << setw(15)
                << ("c" + Foam::name(cloudNo) + varName).c_str() << ' '
                << (masked/ensight::FileName(varName)).c_str() << nl;
        }

        ++cloudNo;
    }


    //
    // TIME
    //

    if (!timesUsed_.empty())
    {
        *os_
            << nl << "TIME" << nl;

        // timeset 1
        const scalar timeCorrection = writeTimeset();

        // timeset geometry
        if (tsGeom > 1)
        {
            writeTimeset(tsGeom, geomTimes_, timeCorrection);
        }

        // timeset cloud
        if (tsCloud > 1)
        {
            writeTimeset(tsCloud, cloudTimes_, timeCorrection);
        }

        *os_
            << "# end" << nl;
    }

    *os_ << flush;
    changed_ = false;
}


Foam::autoPtr<Foam::ensightGeoFile>
Foam::ensightCase::newGeometry
(
    bool moving
) const
{
    autoPtr<ensightGeoFile> filePtr;

    if (UPstream::master())
    {
        // Set the path of the ensight file
        fileName path;

        if (moving)
        {
            // Moving mesh: write as "data/********/geometry"
            path = dataDir()/padded(timeIndex_);
        }
        else
        {
            // Static mesh: write as "data/constant/geometry"
            path = dataDir()/word("constant");
        }
        Foam::mkDir(path);

        noteGeometry(moving);   // note for later use

        filePtr.reset(new ensightGeoFile(path, geometryName, format()));

        // Before 2024-05 also implicitly called beginGeometry()
    }

    return filePtr;
}


Foam::autoPtr<Foam::ensightFile>
Foam::ensightCase::newCloud
(
    const word& cloudName
) const
{
    autoPtr<ensightFile> filePtr;

    if (UPstream::master())
    {
        filePtr = createCloudFile(cloudName, "positions");
        auto& os = filePtr();

        // Tag binary format (just like geometry files)
        os.writeBinaryHeader();

        // Description
        os.write(cloud::prefix/cloudName);
        os.newline();

        noteCloud(cloudName);   // note for later use
    }

    return filePtr;
}


void Foam::ensightCase::rewind() const
{
    if (os_)  // master only
    {
        os_->stdStream().seekp(0, std::ios_base::beg);
    }
}


Foam::Ostream& Foam::ensightCase::printInfo(Ostream& os) const
{
    os  << "Ensight case:" << nl
        << "   path: "   << ensightDir_ << nl
        << "   name: "   << caseName_   << nl
        << "   format: " << format()    << nl;

    if (nodeValues())
    {
        os << "   values per node" << nl;
    }

    return os;
}


// ************************************************************************* //

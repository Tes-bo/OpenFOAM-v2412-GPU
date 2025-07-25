/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "memInfo.H"
#include "IOstreams.H"
#include "OSspecific.H"  // For pid()

#include <cstdlib>
#include <fstream>
#include <string>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::memInfo::memInfo()
:
    peak_(0),
    size_(0),
    rss_(0),
    free_(0)
{
    populate();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::memInfo::good() const noexcept
{
    return peak_ > 0;
}


void Foam::memInfo::clear() noexcept
{
    peak_ = size_ = rss_ = free_ = 0;
}


void Foam::memInfo::populate()
{
    std::string line;

    // "/proc/meminfo"
    // ===========================
    // MemTotal:       65879268 kB
    // MemFree:        51544256 kB
    // MemAvailable:   58999636 kB
    // Buffers:            2116 kB
    // ...
    // Stop parsing when known keys have been extracted
    {
        std::ifstream is("/proc/meminfo");

        for
        (
            unsigned nkeys = 1;
            nkeys && is.good() && std::getline(is, line);
            /*nil*/
        )
        {
            const auto delim = line.find(':');
            if (delim == std::string::npos)
            {
                continue;
            }

            const std::string key(line.substr(0, delim));

            // std::stol() skips whitespace before using as many digits as
            // possible. So just need to skip over the ':' and let stol do
            // the rest

            if (key == "MemFree")
            {
                free_ = std::stol(line.substr(delim+1));
                --nkeys;
            }
        }
    }

    // "/proc/PID/status"
    // ===========================
    // VmPeak:    15920 kB
    // VmSize:    15916 kB
    // VmLck:         0 kB
    // VmPin:         0 kB
    // VmHWM:      6972 kB
    // VmRSS:      6972 kB
    // ...
    // Stop parsing when known keys have been extracted

    // These units are kibi-btyes (1024)
    {
        std::ifstream is("/proc/" + std::to_string(Foam::pid()) + "/status");

        for
        (
            unsigned nkeys = 3;
            nkeys && is.good() && std::getline(is, line);
            /*nil*/
        )
        {
            const auto delim = line.find(':');
            if (delim == std::string::npos)
            {
                continue;
            }

            const std::string key(line.substr(0, delim));

            // std::stoi() skips whitespace before using as many digits as
            // possible. So just need to skip over the ':' and let stoi do
            // the rest

            if (key == "VmPeak")
            {
                peak_ = std::stol(line.substr(delim+1));
                --nkeys;
            }
            else if (key == "VmSize")
            {
                size_ = std::stol(line.substr(delim+1));
                --nkeys;
            }
            else if (key == "VmRSS")
            {
                rss_ = std::stol(line.substr(delim+1));
                --nkeys;
            }
        }
    }
}


const Foam::memInfo& Foam::memInfo::update()
{
    clear();
    populate();
    return *this;
}


void Foam::memInfo::writeEntries(Ostream& os) const
{
    os.writeEntry("size", size_);
    os.writeEntry("peak", peak_);
    os.writeEntry("rss", rss_);
    os.writeEntry("free", free_);
    os.writeEntry("units", "kB");
}


void Foam::memInfo::writeEntry(const word& keyword, Ostream& os) const
{
    os.beginBlock(keyword);
    writeEntries(os);
    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Foam::Istream& Foam::operator>>(Istream& is, memInfo& m)
// {
//     is.readBegin("memInfo");
//     is  >> m.peak_ >> m.size_ >> m.rss_ >> m.free_;
//     is.readEnd("memInfo");
//
//     is.check(FUNCTION_NAME);
//     return is;
// }


Foam::Ostream& Foam::operator<<(Ostream& os, const memInfo& m)
{
    os  << token::BEGIN_LIST
        << m.peak() << token::SPACE
        << m.size() << token::SPACE
        << m.rss()  << token::SPACE
        << m.free()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //

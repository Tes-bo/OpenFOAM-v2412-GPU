/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "cpuInfo.H"
#include "IOstreams.H"

#include <fstream>

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope function
// split things like "a key word\t: value information"
// into ("a_key_word", "value information")

static bool split(const std::string& line, std::string& key, std::string& val)
{
    key.clear();
    val.clear();

    const auto keyLen = line.find_first_of("\t:");
    const auto sep = line.find(':');

    if (keyLen == std::string::npos || sep == std::string::npos)
    {
        return false;
    }

    const auto begVal = line.find_first_not_of(" :", sep);

    if (begVal == std::string::npos)
    {
        return false;
    }

    key = line.substr(0, keyLen);
    val = line.substr(begVal);

    // Avoid spaces in key - replace with '_'
    for (auto iter = key.begin(); iter < key.end(); ++iter)
    {
        if (*iter == ' ')
        {
            *iter = '_';
        }
    }

    // std::cerr<<"key=<" << key << "> val=<" << val << ">\n";

    return true;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Parse the following type of content.
// A TAB separates the keyword from content. Eg,
//
// "cpu cores\t: 6"
//
// ===========================
// processor       : 0
// vendor_id       : GenuineIntel
// cpu family      : 6
// model           : 63
// model name      : Intel(R) Xeon(R) CPU E5-2620 v3 @ 2.40GHz
// stepping        : 2
// microcode       : 0x35
// cpu MHz         : 1200.000
// cache size      : 15360 KB
// physical id     : 0
// siblings        : 12
// core id         : 0
// cpu cores       : 6
// apicid          : 0
// initial apicid  : 0
// fpu             : yes
// fpu_exception   : yes
// cpuid level     : 15
// wp              : yes
// flags           : fpu vme ...
// bugs            :
// bogomips        : 4789.15
// clflush size    : 64
// cache_alignment : 64
// address sizes   : 46 bits physical, 48 bits virtual
// power management:

void Foam::cpuInfo::populate()
{
    int ncpu = 0;
    std::string line, key, val;

    std::ifstream is("/proc/cpuinfo");
    while (is.good() && std::getline(is, line))
    {
        if (!split(line, key, val))
        {
            continue;
        }

        if (key == "processor")
        {
            if (ncpu++)
            {
                break; // stop after the first cpu
            }
        }
        else if (key == "vendor_id")   { vendor_id  = val; }
        else if (key == "model_name")  { model_name = val; }
        else if (key == "cpu_family")  { cpu_family = std::stoi(val); }
        else if (key == "model")       { model = std::stoi(val); }
        else if (key == "cpu_MHz")     { cpu_MHz = std::stof(val); }
        else if (key == "cpu_cores")   { cpu_cores = std::stoi(val);  }
        else if (key == "siblings")    { siblings = std::stoi(val); }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cpuInfo::cpuInfo()
:
    vendor_id(),
    model_name(),
    cpu_family(-1),
    model(-1),
    cpu_MHz(0),
    siblings(0),
    cpu_cores(0)
{
    populate();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cpuInfo::writeEntries(Ostream& os) const
{
    if (!vendor_id.empty())
    {
        os.writeEntry("vendor_id", vendor_id);
    }
    if (!model_name.empty())
    {
        os.writeEntry("model_name", model_name);
    }

    os.writeEntryIfDifferent<int>("cpu_family", -1, cpu_family);
    os.writeEntryIfDifferent<int>("model", -1, model);
    os.writeEntryIfDifferent<float>("cpu_MHz", 0, cpu_MHz);
    os.writeEntryIfDifferent<int>("cpu_cores", 0, cpu_cores);
    os.writeEntryIfDifferent<int>("siblings", 0, siblings);
}


void Foam::cpuInfo::writeEntry(const word& keyword, Ostream& os) const
{
    os.beginBlock(keyword);
    writeEntries(os);
    os.endBlock();
}


// ************************************************************************* //

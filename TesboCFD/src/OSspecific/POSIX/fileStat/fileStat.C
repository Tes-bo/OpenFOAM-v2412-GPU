/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "fileStat.H"
#include "IOstreams.H"
#include "timer.H"

#include <unistd.h>
#ifndef __APPLE__
    #include <sys/sysmacros.h>
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileStat::fileStat()
:
    valid_(false)
{}


Foam::fileStat::fileStat
(
    const char* fName,
    const bool followLink,
    const unsigned int maxTime
)
:
    valid_(false)
{
    if (!fName || !fName[0])
    {
        return;
    }

    // Work on volatile
    volatile bool locIsValid = false;

    timer myTimer(maxTime);

    if (!timedOut(myTimer))
    {
        if (followLink)
        {
            locIsValid = (::stat(fName, &status_) == 0);
        }
        else
        {
            locIsValid = (::lstat(fName, &status_) == 0);
        }
    }

    // Copy into (non-volatile, possible register based) member var
    valid_ = locIsValid;
}


Foam::fileStat::fileStat
(
    const fileName& fName,
    const bool followLink,
    const unsigned int maxTime
)
:
    fileStat(fName.c_str(), followLink, maxTime)
{}


Foam::fileStat::fileStat(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::fileStat::size() const
{
    return valid_ ? label(status_.st_size) : 0;
}


time_t Foam::fileStat::modTime() const
{
    return valid_ ? status_.st_mtime : 0;
}


double Foam::fileStat::dmodTime() const
{
    return
    (
        valid_
      ?
        #ifdef __APPLE__
        (status_.st_mtime + 1e-9*status_.st_mtimespec.tv_nsec)
        #else
        (status_.st_mtime + 1e-9*status_.st_mtim.tv_nsec)
        #endif
      : 0
    );
}


bool Foam::fileStat::sameDevice(const fileStat& other) const
{
    return
        valid_
     && (
            major(status_.st_dev) == major(other.status_.st_dev)
         && minor(status_.st_dev) == minor(other.status_.st_dev)
        );
}


bool Foam::fileStat::sameINode(const fileStat& other) const
{
    return valid_ && (status_.st_ino == other.status_.st_ino);
}


bool Foam::fileStat::sameINode(const label iNode) const
{
    return valid_ && (status_.st_ino == ino_t(iNode));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, fileStat& fs)
{
    FixedList<label, 13> list(is);

    fs.valid_ = list[0];

    dev_t st_dev = makedev(list[1], list[2]);
    fs.status_.st_dev = st_dev;

    fs.status_.st_ino = list[3];
    fs.status_.st_mode = list[4];
    fs.status_.st_uid = list[5];
    fs.status_.st_gid = list[6];

    dev_t st_rdev = makedev(list[7], list[8]);
    fs.status_.st_rdev = st_rdev;

    fs.status_.st_size = list[9];
    fs.status_.st_atime = list[10];
    fs.status_.st_mtime = list[11];
    fs.status_.st_ctime = list[12];

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const fileStat& fs)
{
    FixedList<label, 13> list;

    list[0] = label(fs.valid_);
    list[1] = label(major(fs.status_.st_dev));
    list[2] = label(minor(fs.status_.st_dev));
    list[3] = label(fs.status_.st_ino);
    list[4] = label(fs.status_.st_mode);
    list[5] = label(fs.status_.st_uid);
    list[6] = label(fs.status_.st_gid);
    list[7] = label(major(fs.status_.st_rdev));
    list[8] = label(minor(fs.status_.st_rdev));
    list[9] = label(fs.status_.st_size);
    list[10] = label(fs.status_.st_atime);
    list[11] = label(fs.status_.st_mtime);
    list[12] = label(fs.status_.st_ctime);

    return os << list;
}


// ************************************************************************* //

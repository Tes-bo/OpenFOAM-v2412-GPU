/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "CStringList.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

int Foam::CStringList::reset
(
    std::initializer_list<const char* const> input
)
{
    clear();

    if (!input.size())
    {
        // Special handling of an empty list
        argv_ = new char*[1];
        argv_[0] = nullptr;     // Final nullptr terminator
        return 0;
    }

    // Count overall required string length, including each trailing nul char
    for (const char* const s : input)
    {
        // nbytes_ += Foam::string::length(s) + 1
        nbytes_ += (s ? strlen(s) : 0) + 1;
    }
    --nbytes_;  // Do not include final nul char in overall count

    argv_ = new char*[input.size()+1];  // Extra +1 for terminating nullptr
    data_ = new char[nbytes_+1];        // Extra +1 for terminating nul char

    argv_[0] = data_;   // Starts here

    for (const char* const s : input)
    {
        char *next = stringCopy(argv_[argc_], s);
        argv_[++argc_] = next;  // The start of next string
    }

    argv_[argc_] = nullptr;     // Final nullptr terminator

    return argc_;
}



// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CStringList& list)
{
    const int n = list.size();

    bool separator = false;

    for (int i = 0; i < n; ++i)
    {
        const char* p = list.get(i);

        if (p && *p)
        {
            if (separator) os << ' ';
            os << p;

            separator = true;
        }
    }

    return os;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "genericRagelLemonDriver.H"
#include "error.H"
#include "word.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parsing::genericRagelLemonDriver::genericRagelLemonDriver()
:
    content_(std::cref<std::string>(Foam::string::null)),
    start_(0),
    length_(0),
    position_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::parsing::genericRagelLemonDriver::clear()
{
    content_ = std::cref<std::string>(Foam::string::null);
    start_ = 0;
    length_ = 0;
    position_ = 0;
}


void Foam::parsing::genericRagelLemonDriver::content
(
    const std::string& s,
    size_t pos,
    size_t len
)
{
    content_ = std::cref<std::string>(s);
    start_ = pos;
    length_ = len;
    position_ = 0;
}


std::string::const_iterator
Foam::parsing::genericRagelLemonDriver::cbegin() const
{
    const std::string& s = content_.get();

    if (start_ >= s.length())
    {
        return s.cend();
    }

    return s.cbegin() + start_;
}


std::string::const_iterator
Foam::parsing::genericRagelLemonDriver::cend() const
{
    const std::string& s = content_.get();

    if (length_ == std::string::npos || start_ >= s.length())
    {
        return s.cend();
    }

    const size_t strEnd = start_ + length_;

    if (strEnd >= s.length())
    {
        return s.cend();
    }

    return s.cbegin() + strEnd;
}


Foam::Ostream& Foam::parsing::genericRagelLemonDriver::printBuffer
(
    Ostream& os
) const
{
    const auto endIter = cend();

    for (auto iter = cbegin(); iter != endIter; ++iter)
    {
        char c(*iter);

        // if (!c) break;

        if (c == '\t')
        {
            // Flatten tab to single space for better alignment
            os  << ' ';
        }
        else
        {
            os  << c;
        }
    }

    return os;
}


void Foam::parsing::genericRagelLemonDriver::reportFatal
(
    const std::string& msg
) const
{
    if (position_)
    {
        reportFatal(msg, position_);
    }
    else
    {
        auto& os = FatalIOError
        (
            FUNCTION_NAME,
            __FILE__,
            __LINE__,
            ""
        );

        os  << nl << msg.c_str() << " in expression\n"
            << "<<<<\n";

        printBuffer(os)
            << "\n>>>>\n"
            << exit(Foam::FatalIOError);
    }
}


void Foam::parsing::genericRagelLemonDriver::reportFatal
(
    const std::string& msg,
    size_t pos
) const
{
    auto& os = Foam::FatalIOError
    (
        FUNCTION_NAME,
        __FILE__,
         __LINE__,
        ""
    );

    os  << nl << msg.c_str()
        << " in expression at position:" << long(pos) << nl
        << "<<<<\n";

    const auto begIter = cbegin();
    const auto endIter = cend();

    // Position of newline(s)
    size_t newline0 = 0, newline1 = 0;

    auto iter = begIter;

    for (/*nil*/; iter != endIter; ++iter)
    {
        char c(*iter);

        if ('\t' == c)
        {
            // Flatten tab to single space for better alignment
            os  << ' ';
        }
        else if ('\n' == c)
        {
            os  << c;

            newline1 = (iter-begIter);

            if (newline1 < pos)
            {
                newline0 = newline1;
            }
            else
            {
                ++iter;
                break;
            }
        }
        else
        {
            os  << c;
        }
    }

    if (newline0 == newline1 || newline1 == pos)
    {
        os  << '\n';
    }

    size_t col = std::min(newline0, newline1);
    if (col < pos)
    {
        // This still isn't quite right
        col = pos - col;
        if (col) --col;

        for (/*nil*/; col; --col)
        {
            os  << ' ';
        }
    }

    os  << "^^^^ near here\n";

    // Finish output
    for (/*nil*/; iter != endIter; ++iter)
    {
        char c(*iter);

        if ('\t' == c)
        {
            // Flatten tab to single space for better alignment
            os  << ' ';
        }
        else
        {
            os  << c;
        }
    }

    os  << "\n>>>>\n"
        << exit(Foam::FatalIOError);
}


// ************************************************************************* //

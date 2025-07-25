/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

inline std::string::size_type Foam::string::find_ext(const std::string& str)
{
    const auto i = str.find_last_of("./");

    if (i == npos || i == 0 || str[i] == '/')
    {
        return npos;
    }

    return i;
}


inline std::string::size_type Foam::string::find_ext() const
{
    return find_ext(*this);
}


inline bool Foam::string::has_ext() const
{
    return (std::string::npos != find_ext());
}


inline bool Foam::string::has_ext(const char* ending) const
{
    const auto n1 = size();
    const auto n2 = strlen(ending);

    // Like ends_with with extra check for dot
    return
    (
        n1 > n2
     && operator[](n1-n2-1) == '.'  // Require a dot separator
     && !compare(n1-n2, npos, ending, n2)
    );
}


inline bool Foam::string::has_ext(const std::string& ending) const
{
    const auto n1 = size();
    const auto n2 = ending.size();

    // Like ends_with with extra check for dot
    return
    (
        n1 > n2
     && operator[](n1-n2-1) == '.'  // Require a dot separator
     && !compare(n1-n2, npos, ending)
    );
}


inline bool Foam::string::remove_path()
{
    const auto i = rfind('/');

    if (npos == i)
    {
        return false;
    }

    erase(0, i+1);
    return true;
}


inline bool Foam::string::remove_ext()
{
    const auto i = find_ext();

    if (npos == i)
    {
        return false;
    }

    erase(i);
    return true;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::string::string(const std::string& str)
:
    std::string(str)
{}


inline Foam::string::string(std::string&& str)
:
    std::string(std::move(str))
{}


inline Foam::string::string(const char* str)
:
    std::string(str)
{}


inline Foam::string::string(const char* str, const size_type len)
:
    std::string(str, len)
{}


inline Foam::string::string(const char c)
:
    std::string(1, c)
{}


inline Foam::string::string(const size_type len, const char c)
:
    std::string(len, c)
{}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class StringType>
inline bool Foam::string::valid(const std::string& str)
{
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        if (!StringType::valid(*iter))
        {
            return false;
        }
    }

    return true;
}


template<class StringType>
inline bool Foam::string::stripInvalid(std::string& str)
{
    if (!string::valid<StringType>(str))
    {
        size_type nChar = 0;
        iterator outIter = str.begin();

        for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
        {
            const char c = *iter;

            if (StringType::valid(c))
            {
                *outIter = c;
                ++outIter;
                ++nChar;
            }
        }

        str.erase(nChar);

        return true;
    }

    return false;
}


template<class StringType>
inline StringType Foam::string::validate(const std::string& str)
{
    StringType out;
    out.resize(str.size());

    size_type len = 0;
    for (auto iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        const char c = *iter;
        if (StringType::valid(c))
        {
            out[len] = c;
            ++len;
        }
    }

    out.erase(len);

    return out;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::string::match(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


inline void Foam::string::swap(std::string& str)
{
    if (this != &str)
    {
        // Self-swap is a no-op
        std::string::swap(str);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::string::operator()(const std::string& text) const
{
    return !compare(text);  // Always compare as literal string
}


// ************************************************************************* //

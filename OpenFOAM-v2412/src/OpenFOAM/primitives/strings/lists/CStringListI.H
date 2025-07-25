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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline char* Foam::CStringList::stringCopy(char* dest, const char* src)
{
    dest = std::copy_n(src, (src ? strlen(src) : 0), dest);
    *(dest++) = '\0';
    return dest;
}


inline char* Foam::CStringList::stringCopy(char *dest, const std::string& src)
{
    dest = std::copy_n(src.data(), src.size(), dest);
    *(dest++) = '\0';
    return dest;
}


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline int Foam::CStringList::count(const char * const argv[])
{
    int n = 0;
    if (argv)
    {
        while (argv[n])
        {
            ++n;
        }
    }
    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::CStringList::CStringList() noexcept
:
    argc_(0),
    nbytes_(0),
    argv_(nullptr),
    data_(nullptr)
{}


inline Foam::CStringList::CStringList
(
    std::initializer_list<const char* const> input
)
:
    CStringList()
{
    reset(input);
}


template<class StringType>
inline Foam::CStringList::CStringList(const UList<StringType>& input)
:
    CStringList()
{
    reset(input);
}


template<class StringType>
inline Foam::CStringList::CStringList(const SubStrings<StringType>& input)
:
    CStringList()
{
    reset(input);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::CStringList::~CStringList()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::CStringList::clear()
{
    argc_ = 0;
    nbytes_ = 0;

    if (data_)
    {
        delete[] data_;
        data_ = nullptr;
    }
    if (argv_)
    {
        delete[] argv_;
        argv_ = nullptr;
    }
}


inline bool Foam::CStringList::empty() const noexcept
{
    return !argc_;
}


inline int Foam::CStringList::size() const noexcept
{
    return argc_;
}


inline char** Foam::CStringList::strings() const noexcept
{
    return argv_;
}


inline char** Foam::CStringList::strings(int start) const
{
    return &(argv_[start]);
}


template<class StringType>
inline int Foam::CStringList::reset(const UList<StringType>& input)
{
    return resetContent(input);
}


template<class StringType>
inline int Foam::CStringList::reset(const SubStrings<StringType>& input)
{
    return resetContent(input);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const char* Foam::CStringList::get(int i) const
{
    return argv_[i];
}


inline const char* Foam::CStringList::operator[](int i) const
{
    return argv_[i];
}


// ************************************************************************* //

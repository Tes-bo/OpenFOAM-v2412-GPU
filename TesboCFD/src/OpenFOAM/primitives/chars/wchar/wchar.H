/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Primitive
    wchar_t

Description
    A wide-character and a pointer to a wide-character string.

SourceFiles
    wchar.C

See also
    http://en.wikipedia.org/wiki/UTF-8
    http://en.wikibooks.org/wiki/Unicode/Character_reference

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_wchar_H
#define Foam_primitives_wchar_H

#include <cwchar>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Output wide character (Unicode) as UTF-8
Ostream& operator<<(Ostream& os, const wchar_t wc);

//- Output wide character (Unicode) string as UTF-8
Ostream& operator<<(Ostream& os, const wchar_t* wstr);

//- Output wide character (Unicode) string as UTF-8
Ostream& operator<<(Ostream& os, const std::wstring& wstr);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

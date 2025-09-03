/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

Class
    Foam::functionEntries::inputMode

Description
    Specify the global input mode when reading further dictionaries,
    expects a single word to follow.

    An example of \c \#inputMode directive:
    \verbatim
        #inputMode merge
    \endverbatim

    The possible input modes:
      - \par merge
        merge sub-dictionaries when possible
      - \par overwrite
        keep last entry and silently remove previous ones
      - \par protect
        keep initial entry and silently ignore subsequent ones
      - \par warn
        keep initial entry and warn about subsequent ones
      - \par error
        issue a FatalError for duplicate entries
      - \par default
        The default treatment - currently identical to \b merge.

SourceFiles
    inputMode.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_inputMode_H
#define functionEntries_inputMode_H

#include "entry.H"
#include "functionEntry.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                          Class inputMode Declaration
\*---------------------------------------------------------------------------*/

class inputMode
:
    public functionEntry
{
public:

    //- Text representations of the selectable input modes
    static const Enum<entry::inputMode> selectableNames;

    //- Change the global input-mode.
    static bool execute(dictionary& unused, Istream& is);
};


/*---------------------------------------------------------------------------*\
                      Class inputModeDefault Declaration
\*---------------------------------------------------------------------------*/

//- Temporarily change inputMode to %protect for the following entry
class inputModeDefault
:
    public functionEntry
{
public:

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


/*---------------------------------------------------------------------------*\
                       Class inputModeMerge Declaration
\*---------------------------------------------------------------------------*/

//- Temporarily change inputMode to %merge for the following entry
class inputModeMerge
:
    public functionEntry
{
public:

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


/*---------------------------------------------------------------------------*\
                     Class inputModeOverwrite Declaration
\*---------------------------------------------------------------------------*/

//- Temporarily change inputMode to %overwrite for the following entry
class inputModeOverwrite
:
    public functionEntry
{
public:

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


/*---------------------------------------------------------------------------*\
                        Class inputModeWarn Declaration
\*---------------------------------------------------------------------------*/

//- Temporarily change inputMode to %warn for the following entry
class inputModeWarn
:
    public functionEntry
{
public:

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


/*---------------------------------------------------------------------------*\
                       Class inputModeError Declaration
\*---------------------------------------------------------------------------*/

//- Temporarily change inputMode to %error for the following entry
class inputModeError
:
    public functionEntry
{
public:

    //- Execute in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

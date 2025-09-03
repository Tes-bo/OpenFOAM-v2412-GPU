/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::functionEntries::includeFuncEntry

Description
    Specify a functionObject dictionary file to include, expects the
    functionObject name to follow with option arguments (without quotes).

    Uses functionObjectList::readFunctionObject() method to search
    within user/group/other "caseDicts/postProcessing" directories
    as per Foam::functionObjectList::findDict() and Foam::findEtcFiles()

    The optional field arguments included in the name are inserted in 'field' or
    'fields' entries in the functionObject dictionary and included in the name
    of the functionObject entry to avoid conflict.

    Examples:
    \verbatim
        #includeFunc Q
        #includeFunc components(U)
        #includeFunc mag(Ux)
        #includeFunc mag(p)
    \endverbatim

    Other dictionary entries may also be specified using named arguments, for
    example the \c name of the \c faceZone in the \c flowRateFaceZone \c
    functionObject configuration set and the \c orientedFields entry which
    defaults to \c phi may also be overridden as required, e.g.

    \verbatim
        #includeFunc flowRateFaceZone(name=fZone1)
        #includeFunc flowRateFaceZone(orientedFields=(phiAlpha),name=fZone1)
    \endverbatim

See also
    Foam::functionObjectList

SourceFiles
    includeFuncEntry.C

\*---------------------------------------------------------------------------*/

#ifndef functionEntries_includeFuncEntry_H
#define functionEntries_includeFuncEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeFuncEntry Declaration
\*---------------------------------------------------------------------------*/

class includeFuncEntry
:
    public functionEntry
{
public:

    //- Run function in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fileFormats::AC3DsurfaceFormatCore

Description
    Internal class used by the AC3DsurfaceFormat

SourceFiles
    AC3DsurfaceFormatCore.C

\*---------------------------------------------------------------------------*/

#ifndef AC3DsurfaceFormatCore_H
#define AC3DsurfaceFormatCore_H

#include "Fstream.H"
#include "surfZone.H"
#include "StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
             Class fileFormats::AC3DsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class AC3DsurfaceFormatCore
{
protected:

    // Protected Static Member Functions

        //- Read a type from string
        template<class Type>
        static Type parse(const string& str)
        {
            IStringStream is(str);

            Type t;
            is >> t;

            return t;
        }

        //- Read cmd, args from IFstream
        //  The cmd is the content up to the first space, args is the balance
        //  of the line (after the space).
        //  \return if the command was read
        static bool readCmd(IFstream& is, string& cmd, string& args);

        //- Read up to a line starting with cmd. Sets args to rest of line.
        //
        //  \return true if found, false if stream is not good anymore.
        static bool cueTo(IFstream&, const string& cmd, string& args);

        //- Like cueTo(), but FatalError if not found.
        //  \return the command args on success
        static string cueToOrDie
        (
            IFstream&,
            const string& cmd,
            const string& errorMsg=string::null
        );

        //- Write header with materials for each zone
        static void writeHeader(Ostream& os, const UList<surfZone>& zones);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

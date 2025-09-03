/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "blockMeshTools.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static inline const Foam::entry* resolveLabel(const entry& e, const label val)
{
    return
    (
        (e.isStream() && e.stream().front().isLabel(val))
      ? &e : nullptr
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockMeshTools::read
(
    Istream& is,
    label& val,
    const dictionary& dict
)
{
    token t(is);
    if (t.isLabel())
    {
        val = t.labelToken();
    }
    else if (t.isWord())
    {
        const word& varName = t.wordToken();
        const entry* eptr =
            dict.findScoped(varName, keyType::REGEX_RECURSIVE);

        if (eptr)
        {
            // Read as label
            val = Foam::readLabel(eptr->stream());
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Undefined variable "
                << varName << ". Valid variables are " << dict
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "Illegal token " << t.info()
            << " when trying to read label"
            << exit(FatalIOError);
    }

    is.fatalCheck(FUNCTION_NAME);
}


Foam::label Foam::blockMeshTools::read
(
    Istream& is,
    const dictionary& dict
)
{
    label val(0);
    read(is, val, dict);
    return val;
}


void Foam::blockMeshTools::write
(
    Ostream& os,
    const label val,
    const dictionary& dict
)
{
    for (const entry& e : dict)
    {
        const entry* eptr = resolveLabel(e, val);
        if (eptr)
        {
            os << eptr->keyword();
            return;
        }
    }
    os << val;
}


const Foam::entry* Foam::blockMeshTools::findEntry
(
    const dictionary& dict,
    const label val
)
{
    for (const entry& e : dict)
    {
        const entry* eptr = resolveLabel(e, val);
        if (eptr)
        {
            return eptr;
        }
    }

    return nullptr;
}


// ************************************************************************* //

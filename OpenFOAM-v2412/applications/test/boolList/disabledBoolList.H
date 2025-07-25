/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::disabledBoolList

Description
    A trivial structure to use as a boolList replacement when testing
    compilation without using the [] accessors.

    Example,
    \code
    const disableBoolList blockedFace;
    ...

    if (blockedFace.test(facei)) ...    // Good
    if (blockedFace[facei]) ...         // Compile error
    \endcode

\*---------------------------------------------------------------------------*/

#ifndef Foam_disabledBoolList_H
#define Foam_disabledBoolList_H

/*---------------------------------------------------------------------------*\
                      Class disabledBoolList Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{

struct disabledBoolList
{
    bool empty() const { return true; }

    int size() const { return 0; }

    bool test(int) const { return true; }

    void set(bool) {}

    // Perhaps not?
    /// bool operator()(const label i) const { return true; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

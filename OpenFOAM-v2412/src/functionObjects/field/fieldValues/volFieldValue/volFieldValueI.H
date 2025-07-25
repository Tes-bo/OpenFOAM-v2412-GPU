/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

inline bool Foam::functionObjects::fieldValues::volFieldValue::is_magOp()
const noexcept
{
    // Operation tagged to use mag
    return (operation_ & typeAbsolute);
}


inline bool Foam::functionObjects::fieldValues::volFieldValue::is_weightedOp()
const noexcept
{
    // Operation tagged to require a weight field
    return (operation_ & typeWeighted);
}


inline bool Foam::functionObjects::fieldValues::volFieldValue::canWeight
(
    const scalarField& fld
) const
{
    // Non-empty on some processor
    return returnReduceOr(!fld.empty());
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::externalFileCoupler::initialized() const
{
    return runState_ != NONE;
}


inline bool Foam::externalFileCoupler::slaveFirst() const
{
    return slaveFirst_;
}


inline const Foam::fileName& Foam::externalFileCoupler::commDirectory() const
{
    return commsDir_;
}


inline Foam::fileName Foam::externalFileCoupler::resolveFile
(
    const word& file
) const
{
    return fileName(commDirectory()/file);
}


inline Foam::fileName Foam::externalFileCoupler::lockFile() const
{
    return resolveFile(lockName + ".lock");
}


// ************************************************************************* //

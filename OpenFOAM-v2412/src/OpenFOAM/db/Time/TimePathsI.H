/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::TimePaths::TimePaths
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemDirName,
    const word& constantDirName
)
:
    TimePaths
    (
        false,              // processorCase
        rootPath,           // rootPath
        false,              // distributed
        caseName,           // globalCaseName
        caseName,           // caseName
        systemDirName,
        constantDirName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::TimePaths::distributed() const noexcept
{
    return distributed_;
}


inline bool Foam::TimePaths::processorCase() const noexcept
{
    return processorCase_;
}


inline bool Foam::TimePaths::processorCase(bool on) noexcept
{
    bool old(processorCase_);
    processorCase_ = on;
    return old;
}


inline const Foam::fileName& Foam::TimePaths::rootPath() const noexcept
{
    return rootPath_;
}


inline const Foam::fileName& Foam::TimePaths::globalCaseName() const noexcept
{
    return globalCaseName_;
}


inline const Foam::fileName& Foam::TimePaths::caseName() const noexcept
{
    return case_;
}


inline Foam::fileName& Foam::TimePaths::caseName() noexcept
{
    return case_;
}


inline Foam::fileName Foam::TimePaths::path() const
{
    return rootPath()/caseName();
}


inline Foam::fileName Foam::TimePaths::globalPath() const
{
    return rootPath()/globalCaseName();
}


inline Foam::fileName Foam::TimePaths::relativePath
(
    const fileName& input,
    const bool caseTag
) const
{
    return input.relative(globalPath(), caseTag);
}


inline const Foam::word& Foam::TimePaths::constant() const noexcept
{
    return constant_;
}


inline const Foam::word& Foam::TimePaths::system() const noexcept
{
    return system_;
}


inline Foam::fileName Foam::TimePaths::caseConstant() const
{
    if (processorCase_)
    {
        return ".."/constant();
    }

    return constant();
}


inline Foam::fileName Foam::TimePaths::caseSystem() const
{
    if (processorCase_)
    {
        return ".."/system();
    }

    return system();
}


inline Foam::fileName Foam::TimePaths::constantPath() const
{
    return path()/constant();
}


inline Foam::fileName Foam::TimePaths::systemPath() const
{
    return path()/system();
}


// ************************************************************************* //

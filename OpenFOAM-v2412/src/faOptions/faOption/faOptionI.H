/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
------------------------------------------------------------------------------
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

inline const Foam::word& Foam::fa::option::name() const noexcept
{
    return name_;
}


inline const Foam::fvMesh& Foam::fa::option::mesh() const noexcept
{
    return mesh_;
}


inline const Foam::dictionary& Foam::fa::option::coeffs() const noexcept
{
    return coeffs_;
}


inline bool Foam::fa::option::active() const noexcept
{
    return active_;
}


inline bool Foam::fa::option::active(const bool on) noexcept
{
    bool old(active_);
    active_ = on;
    return old;
}


inline void Foam::fa::option::setApplied(const label fieldi)
{
    applied_[fieldi] = true;
}


inline const Foam::word& Foam::fa::option::regionName() const noexcept
{
    return regionName_;
}


inline const Foam::faMesh& Foam::fa::option::regionMesh() const
{
    if (!regionMeshPtr_)
    {
        regionMeshPtr_.reset(new faMesh(mesh_));
    }
    return *regionMeshPtr_;
}


inline const Foam::volSurfaceMapping& Foam::fa::option::vsm() const
{
    if (!vsmPtr_)
    {
        vsmPtr_.reset(new volSurfaceMapping(this->regionMesh()));
    }
    return *vsmPtr_;
}


// ************************************************************************* //

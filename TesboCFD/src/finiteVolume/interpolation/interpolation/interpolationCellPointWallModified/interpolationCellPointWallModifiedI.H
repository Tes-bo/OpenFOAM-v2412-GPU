/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

template<class Type>
inline Type Foam::interpolationCellPointWallModified<Type>::interpolate
(
    const cellPointWeight& cpw
) const
{
    return interpolationCellPoint<Type>::interpolate(cpw);
}


template<class Type>
inline Type Foam::interpolationCellPointWallModified<Type>::interpolate
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    return interpolationCellPoint<Type>::interpolate
    (
        cellPointWeightWallModified
        (
            this->pMesh_,
            position,
            celli,
            facei
        )
    );
}


template<class Type>
inline Type Foam::interpolationCellPointWallModified<Type>::interpolate
(
    const barycentric& coordinates,
    const tetIndices& tetIs,
    const label facei
) const
{
    if (facei >= 0)
    {
        if (facei != tetIs.face())
        {
            FatalErrorInFunction
                << "specified face " << facei << " inconsistent with the face "
                << "stored by tetIndices: " << tetIs.face()
                << exit(FatalError);
        }

        if (cellPointWeightWallModified::onWall(this->pMesh_, facei))
        {
            return this->psi_[tetIs.cell()];
        }
    }

    // If the wall face selection did not return, then use the normal
    // interpolate method

    return interpolationCellPoint<Type>::interpolate
    (
        coordinates,
        tetIs,
        facei
    );
}


// ************************************************************************* //

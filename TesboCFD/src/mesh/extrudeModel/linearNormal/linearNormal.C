/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "linearNormal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linearNormal, 0);

addToRunTimeSelectionTable(extrudeModel, linearNormal, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearNormal::linearNormal(const dictionary& dict)
:
    extrudeModel(typeName, dict),
    thickness_(coeffDict_.get<scalar>("thickness")),
    firstCellThickness_
    (
        coeffDict_.getOrDefault<scalar>("firstCellThickness", 0)
    ),
    layerPoints_(nLayers_)
{
    if (thickness_ <= 0)
    {
        FatalErrorInFunction
            << "thickness should be positive : " << thickness_
            << exit(FatalError);
    }

    if (nLayers_ > 1 && firstCellThickness_ > 0)
    {
        if (thickness_ <= firstCellThickness_)
        {
            FatalErrorInFunction
                << "firstCellThickness leave no room for further layers"
                << exit(FatalError);
        }

        layerPoints_[0] = firstCellThickness_;

        for (label layer = 1; layer < nLayers_; ++layer)
        {
            layerPoints_[layer] =
                (thickness_ - layerPoints_[0])
                *sumThickness(layer) + layerPoints_[0];
        }
    }
    else
    {
        for (label layer = 0; layer < nLayers_; ++layer)
        {
            layerPoints_[layer] = thickness_*sumThickness(layer + 1);
        }
    }
}


// * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * * //

point linearNormal::operator()
(
    const point& surfacePoint,
    const vector& surfaceNormal,
    const label layer
) const
{
    if (layer <= 0)
    {
        return surfacePoint;
    }

    return surfacePoint + layerPoints_[layer - 1]*surfaceNormal;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// ************************************************************************* //

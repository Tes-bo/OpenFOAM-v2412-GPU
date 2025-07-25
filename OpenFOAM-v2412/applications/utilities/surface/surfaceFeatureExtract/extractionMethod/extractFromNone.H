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

Class
    Foam::surfaceFeaturesExtraction::extractFromNone

Description
    Run-time selectable surface feature extraction - no extraction.
    Primarily useful with self-intersection methods.

    Selectable as "none".
    Optional dictionary entries: "includedAngle", "geometricTestOnly".

SourceFiles
    extractFromNone.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFeaturesExtraction_extractFromNone_H
#define surfaceFeaturesExtraction_extractFromNone_H

#include "surfaceFeaturesExtraction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{

/*---------------------------------------------------------------------------*\
         Class surfaceFeaturesExtraction::extractFromNone Declaration
\*---------------------------------------------------------------------------*/

class extractFromNone
:
    public method
{
public:

    //- Construct from dictionary
    explicit extractFromNone(const dictionary& dict);

    //- Destructor
    virtual ~extractFromNone() = default;

    //- Extracted features from surface (no-op)
    autoPtr<surfaceFeatures> features(const triSurface& surf) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFeaturesExtraction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

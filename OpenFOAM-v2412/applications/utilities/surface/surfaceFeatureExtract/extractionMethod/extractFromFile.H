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
    Foam::surfaceFeaturesExtraction::extractFromFile

Description
    Run-time selectable surface feature extraction.

    Selectable as "extractFromFile".
    Mandatory dictionary entries: "featureEdgeFile".
    Optional dictionary entries: "geometricTestOnly".

SourceFiles
    extractFromFile.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFeaturesExtraction_extractFromFile_H
#define surfaceFeaturesExtraction_extractFromFile_H

#include "surfaceFeaturesExtraction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{

/*---------------------------------------------------------------------------*\
         Class surfaceFeaturesExtraction::extractFromFile Declaration
\*---------------------------------------------------------------------------*/

class extractFromFile
:
    public method
{
    fileName featureEdgeFile_;


public:

    //- Construct from dictionary
    extractFromFile(const dictionary& dict);

    //- Destructor
    virtual ~extractFromFile() = default;

    //- Features loaded (extracted) from featureEdgeFile
    autoPtr<surfaceFeatures> features(const triSurface& surf) const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFeaturesExtraction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

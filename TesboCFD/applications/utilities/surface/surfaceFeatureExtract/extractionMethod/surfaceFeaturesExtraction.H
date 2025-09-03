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

Namespace
    Foam::surfaceFeaturesExtraction

Description
    Namespace for run-time selectable surface feature extraction methods.

Class
    Foam::surfaceFeaturesExtraction::method

Description
    Abstract base for run-time selectable surface feature extraction methods.

SourceFiles
    surfaceFeaturesExtraction.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFeaturesExtraction_method_H
#define surfaceFeaturesExtraction_method_H

#include "surfaceFeatures.H"
#include "dictionary.H"
#include "Switch.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFeaturesExtraction
{

/*---------------------------------------------------------------------------*\
              Class surfaceFeaturesExtraction::method Declaration
\*---------------------------------------------------------------------------*/

class method
{
protected:

    scalar includedAngle_;
    Switch geometricTestOnly_;

    //- Default construct
    method();


public:

    //- Runtime type information
    ClassNameNoDebug("method");


    // Constructors

        //- Construct from dictionary
        explicit method(const dictionary& dict);


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        method,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<method> New(const dictionary& dict);


    //- Destructor
    virtual ~method() = default;


    // Member Functions

        //- The included angle, if set
        scalar includedAngle() const noexcept { return includedAngle_; }

        //- Use geometric test only
        Switch geometricTestOnly() const noexcept { return geometricTestOnly_; }

        //- Extracted features
        virtual autoPtr<surfaceFeatures> features
        (
            const triSurface& surf
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFeaturesExtraction
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

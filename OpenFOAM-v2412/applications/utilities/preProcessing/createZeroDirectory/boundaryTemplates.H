/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::boundaryTemplates

Description
    Class to store boundary template specifications

    Templates are typically stored centrally, and constructed hierarchically.
    The main use is to convert the (user) specified conditions into
    a form that can be inserted into each field file as dictionary entries.

\*---------------------------------------------------------------------------*/

#ifndef boundaryTemplates_H
#define boundaryTemplates_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;

/*---------------------------------------------------------------------------*\
                      Class boundaryTemplates Declaration
\*---------------------------------------------------------------------------*/

class boundaryTemplates
{
    // Private data

        //- Dictionary of boundary templates
        dictionary templates_;

        //- Dictionary of boundary template options
        dictionary options_;


public:

    //- Constructor
    boundaryTemplates
    (
        const fileName& baseDir,
        const Time& runTime,
        const word& solverType
    );


    // Public member functions

        //- Return the dictionary of boundary templates
        const dictionary& templates() const;

        //- Generate a dictionary representation of patch boundary condition
        dictionary generatePatchDict
        (
            const word& regionPrefix,
            const word& fieldName,
            const word& condition,
            const word& category,
            const word& patchType,
            const dictionary& conditionOptions
        ) const;

        //- Check that user supplied patch info is valid
        void checkPatch
        (
            const word& regionPrefix,
            const word& condition,
            const word& category,
            const word& patchType
        ) const;

        //- Return true if condition requires additional user options
        bool optionsRequired
        (
            const word& regionPrefix,
            const word& category,
            const word& patchType
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

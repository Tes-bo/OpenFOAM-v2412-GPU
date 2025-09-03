/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::functionObjects::volFieldSelection

Description
    Helper class to manage solver field selections

SourceFiles
    solverFieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_solverFieldSelection_H
#define Foam_functionObjects_solverFieldSelection_H

#include "fieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class solverFieldSelection Declaration
\*---------------------------------------------------------------------------*/

class solverFieldSelection
:
    public fieldSelection
{
public:

    //- No copy construct
    solverFieldSelection(const solverFieldSelection&) = delete;


    //- Construct from object registry
    explicit solverFieldSelection
    (
        const objectRegistry& obr,
        const bool includeComponents = false
    );


    //- Destructor
    virtual ~solverFieldSelection() = default;


    // Member Functions

        //- Update the selection using current contents of obr_
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

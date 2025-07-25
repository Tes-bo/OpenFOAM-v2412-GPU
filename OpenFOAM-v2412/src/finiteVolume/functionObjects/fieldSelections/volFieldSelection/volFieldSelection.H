/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Helper class to manage volume field selections

SourceFiles
    volFieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_volFieldSelection_H
#define Foam_functionObjects_volFieldSelection_H

#include "fieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class volFieldSelection Declaration
\*---------------------------------------------------------------------------*/

class volFieldSelection
:
    public fieldSelection
{
protected:

    // Protected Member Functions

        //- Add registered GeometricField types to selection
        template<template<class> class PatchType, class MeshType>
        void addRegisteredGeoFields(DynamicList<fieldInfo>& set) const;


public:

    //- No copy construct
    volFieldSelection(const volFieldSelection&) = delete;


    //- Construct from object registry
    explicit volFieldSelection
    (
        const objectRegistry& obr,
        const bool includeComponents = false
    );


    //- Destructor
    virtual ~volFieldSelection() = default;


    // Member Functions

        //- Update the selection using current contents of obr_
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "volFieldSelectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

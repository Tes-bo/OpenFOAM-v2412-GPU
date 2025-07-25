/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::functionObjects::fileFieldSelection

Description
    Helper class to manage file-based field selections

SourceFiles
    fieldSelection.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_fileFieldSelection_H
#define Foam_functionObjects_fileFieldSelection_H

#include "fieldSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class IOobjectList;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class fileFieldSelection Declaration
\*---------------------------------------------------------------------------*/

class fileFieldSelection
:
    public fieldSelection
{
    // Private Member Functions

        //- Add objects of a given type
        template<class Type>
        void addFromFile
        (
            const IOobjectList& objects,
            DynamicList<fieldInfo>& set
        ) const;

        //- Add registered GeometricField types to selection
        template<template<class> class PatchType, class MeshType>
        void addGeoFieldTypes
        (
            const IOobjectList& objects,
            DynamicList<fieldInfo>& set
        ) const;

        //- Add registered Internal types to selection
        void addInternalFieldTypes
        (
            const IOobjectList& objects,
            DynamicList<fieldInfo>& set
        ) const;

        //- Add registered uniform types to selection
        void addUniformFieldTypes
        (
            const IOobjectList& objects,
            DynamicList<fieldInfo>& set
        ) const;


public:

    //- No copy construct
    fileFieldSelection(const fileFieldSelection&) = delete;


    // Constructors

        //- Construct from object registry
        explicit fileFieldSelection
        (
            const objectRegistry& obr,
            const bool includeComponents = false
        );


    //- Destructor
    virtual ~fileFieldSelection() = default;


    // Member Functions

        //- Update the selection
        virtual bool updateSelection();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

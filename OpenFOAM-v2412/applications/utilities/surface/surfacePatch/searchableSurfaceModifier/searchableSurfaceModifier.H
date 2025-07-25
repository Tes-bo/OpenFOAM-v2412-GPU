/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::searchableSurfaceModifier

Description
    Changing a surface

SourceFiles
    searchableSurfaceModifier.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceModifier_H
#define searchableSurfaceModifier_H

#include "dictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;
class searchableSurface;
class searchableSurfaces;

/*---------------------------------------------------------------------------*\
                  Class searchableSurfaceModifier Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceModifier
{
protected:

    // Protected data

        const searchableSurfaces& geometry_;

        //- Input dictionary
        const dictionary dict_;

public:

    //- Runtime type information
    TypeName("searchableSurfaceModifier");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            searchableSurfaceModifier,
            dictionary,
            (
                const searchableSurfaces& geometry,
                const dictionary& dict
            ),
            (geometry, dict)
        );


    // Constructors

        //- Construct from dictionary
        searchableSurfaceModifier(const searchableSurfaces&, const dictionary&);

        //- Clone
        autoPtr<searchableSurfaceModifier> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a reference to the selected searchableSurfaceModifier
        static autoPtr<searchableSurfaceModifier> New
        (
            const word& type,
            const searchableSurfaces&,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceModifier();


    // Member Functions

        //- Do any operation on surface. Return true if anything changed
        virtual bool modify(const labelList&, searchableSurface&) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

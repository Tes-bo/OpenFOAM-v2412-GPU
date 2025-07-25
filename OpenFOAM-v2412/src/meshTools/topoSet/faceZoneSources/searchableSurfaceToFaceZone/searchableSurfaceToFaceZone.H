/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::searchableSurfaceToFaceZone

Description
    A \c topoSetSource to select all faces whose cell-cell centre vector
    intersects with a given \c searchableSurface.

    Operands:
    \table
      Operand   | Type       | Location
      input     | triSurface | $FOAM_CASE/constant/triSurface/\<file\>
      output 1  | faceSet    | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 2  | faceZone   | $FOAM_CASE/constant/polyMesh/faceZones
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceZoneSet;
        action      <action>;

        // Mandatory entries
        source      searchableSurfaceToFaceZone;
        surfaceType <surfaceTypeName>;

        // Optional entries
        surfaceName <surfaceName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Req'd | Dflt
      name       | Name of faceZone                     | word |  yes  | -
      type       | Type name: faceZone                  | word |  yes  | -
      action     | Action applied on faces - see below  | word |  yes  | -
      source     | Source name: searchableSurfaceToFaceZone | word | yes | -
      surfaceType | The searchable surface type        | word |  yes | -
      surfaceName | Name for the IOobject              | word | no   | mesh-name
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceZone from selected faces
      add      | Add selected faces into this faceZone
      subtract | Remove selected faces from this faceZone
    \endverbatim

See also
    - Foam::topoSetSource

SourceFiles
    searchableSurfaceToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_searchableSurfaceToFaceZone_H
#define Foam_searchableSurfaceToFaceZone_H

#include "topoSetFaceZoneSource.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class searchableSurfaceToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceToFaceZone
:
    public topoSetFaceZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Surface
        autoPtr<searchableSurface> surfacePtr_;


public:

    //- Runtime type information
    TypeName("searchableSurfaceToFaceZone");


    // Constructors

        //- Construct surface-type from dictionary
        searchableSurfaceToFaceZone
        (
            const word& surfaceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from dictionary
        searchableSurfaceToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceToFaceZone() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

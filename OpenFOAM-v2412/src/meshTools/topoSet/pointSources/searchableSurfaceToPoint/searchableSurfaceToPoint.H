/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::searchableSurfaceToPoint

Description
    A \c topoSetPointSource to select points which are
    centre enclosed by a given \c searchableSurface.

    Operands:
    \table
      Operand   | Type       | Location
      input     | triSurface | $FOAM_CASE/constant/triSurface/\<file\>
      output    | pointSet   | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        pointSet;
        action      <action>;

        // Mandatory entries
        source      searchableSurfaceToPoint;
        surfaceType <surfaceTypeName>;

        // Optional entries
        surfaceName <surfaceName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                           | Type | Req'd | Dflt
      name       | Name of pointSet                      | word |  yes  | -
      type       | Type name: pointSet                   | word |  yes  | -
      action     | Action applied on points - see below  | word |  yes  | -
      source     | Source name: searchableSurfaceToPoint | word |  yes | -
      surfaceType | The searchable surface type          | word |  yes  | -
      surfaceName | Name for the IOobject                | word | no | mesh-name
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointSet from selected points
      add      | Add selected points into this pointSet
      subtract | Remove selected points from this pointSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    searchableSurfaceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_searchableSurfaceToPoint_H
#define Foam_searchableSurfaceToPoint_H

#include "topoSetPointSource.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class searchableSurfaceToPoint Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceToPoint
:
    public topoSetPointSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- The searchableSurface
        autoPtr<searchableSurface> surf_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;

public:

    //- Runtime type information
    TypeName("searchableSurfaceToPoint");


    // Constructors

        //- Construct surface-type from dictionary
        searchableSurfaceToPoint
        (
            const word& surfaceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from dictionary
        searchableSurfaceToPoint
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceToPoint() = default;


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

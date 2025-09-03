/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::shapeToCell

Description
    A \c topoSetCellSource to select cells
    based on the type of their cell shapes.

    Handles all types of cell shapes known from static collection
    in the \c cellModel (as shown below) and \c splitHex with
    10 degrees feature angle.

    Operands:
    \table
      Operand   | Type    | Location
      output    | cellSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        cellSet;
        action      <action>;

        // Mandatory entries
        source      shapeToCell;
        shape       <shapeTypeName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: shapeToCell            | word |  yes  | -
      shape      | Type of cell shape - see below      | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Options for the \c shape entry:
    \verbatim
      splitHex   | Split-hexahedra cells
      hex        | Hexahedra cells
      wedge      | Wedge cells
      tetWedge   | Tetrahedra-wedge cells
      prism      | Prism cells
      pyr        | Pyramid cells
      tet        | Tetrahedra cells
    \endverbatim

Note
    \c splitHex hardcoded with internal angle < 10 degrees.

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    shapeToCell.C

\*---------------------------------------------------------------------------*/

#ifndef shapeToCell_H
#define shapeToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class shapeToCell Declaration
\*---------------------------------------------------------------------------*/

class shapeToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of cell shape/type
        word shape_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("shapeToCell");


    // Static Data

        //- Cos of feature angle for polyHedral to be splitHex
        static scalar featureCos;


    // Constructors

        //- Construct from components
        shapeToCell(const polyMesh& mesh, const word& shapeName);

        //- Construct from dictionary
        shapeToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        shapeToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~shapeToCell() = default;


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

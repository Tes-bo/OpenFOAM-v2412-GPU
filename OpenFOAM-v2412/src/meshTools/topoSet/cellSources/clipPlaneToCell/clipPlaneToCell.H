/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::clipPlaneToCell

Description
    A \c topoSetCellSource to select all cells
    whose cell centre is above the specified plane

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
        source      clipPlaneToCell;

        point  (x y z);
        normal (x y z);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: clipPlaneToCell        | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                        | Type   | Req'd  | Dflt
      point    | The reference point for the plane  | vector | yes    | -
      normal   | The normal of the plane            | vector | yes    | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    clipPlaneToCell.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_clipPlaneToCell_H
#define Foam_clipPlaneToCell_H

#include "topoSetCellSource.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class clipPlaneToCell Declaration
\*---------------------------------------------------------------------------*/

class clipPlaneToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Point on the plane
        const vector point_;

        //- Normal to the plane
        const vector normal_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("clipPlaneToCell");


    // Constructors

        //- No default construct
        clipPlaneToCell() = delete;

        //- Construct from components
        clipPlaneToCell
        (
            const polyMesh& mesh,
            const point& basePoint,
            const vector& normal
        );

        //- Construct from dictionary
        clipPlaneToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from stream - a single box.
        clipPlaneToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~clipPlaneToCell() = default;


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

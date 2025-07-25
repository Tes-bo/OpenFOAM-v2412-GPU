/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::boundaryToCell

Description
    A \c topoSetCellSource to select all cells possessing
    at least a single face on all external boundaries.

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
        source      boundaryToCell;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Req'd | Dflt
      name       | Name of cellSet                      | word |  yes  | -
      type       | Type name: cellSet                   | word |  yes  | -
      action     | Action applied on cells - see below  | word |  yes  | -
      source     | Source name: boundaryToCell          | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource
    - Foam::patchToCell

SourceFiles
    boundaryToCell.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryToCell_H
#define boundaryToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class boundaryToCell Declaration
\*---------------------------------------------------------------------------*/

class boundaryToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("boundaryToCell");


    // Constructors

        //- Construct from components
        boundaryToCell(const polyMesh& mesh);

        //- Construct from dictionary
        boundaryToCell(const polyMesh& mesh, const dictionary& unused);

        //- Construct from Istream
        boundaryToCell(const polyMesh& mesh, Istream& unused);


    //- Destructor
    virtual ~boundaryToCell() = default;


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

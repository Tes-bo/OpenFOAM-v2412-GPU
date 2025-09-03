/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::targetVolumeToCell

Description
    A \c topoSetCellSource to select cells based on a target
    volume of cells. Adapts a plane until it has enough.

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
        source      targetVolumeToCell;
        volume      10;
        normal      (1 0 0);

        // Optional entries
        set         <cellSetName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: targetVolumeToCell     | word |  yes  | -
      volume     | The target volume [m3]              | scalar | yes | -
      normal     | The plane normal                    | vector | yes | -
      set        | Restrict to named cellSet           | word |  no   | ""
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

SourceFiles
    targetVolumeToCell.C

\*---------------------------------------------------------------------------*/

#ifndef targetVolumeToCell_H
#define targetVolumeToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class bitSet;

/*---------------------------------------------------------------------------*\
                         Class targetVolumeToCell Declaration
\*---------------------------------------------------------------------------*/

class targetVolumeToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Target volume
        const scalar vol_;

        //- Normal of plane to sweep
        const vector normal_;

        //- Optional name of cellSet to calculate volume in
        const word maskSetName_;


    // Private Member Functions

        //- The volume of the selected cells
        scalar volumeOfSet(const bitSet& selected) const;

        label selectCells
        (
            const scalar normalComp,
            const bitSet&,
            bitSet& selected
        ) const;

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("targetVolumeToCell");


    // Constructors

        //- Construct from components
        targetVolumeToCell
        (
            const polyMesh& mesh,
            const scalar vol,
            const vector& normal,
            const word& maskSetName = ""
        );

        //- Construct from dictionary
        targetVolumeToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        targetVolumeToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~targetVolumeToCell() = default;


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

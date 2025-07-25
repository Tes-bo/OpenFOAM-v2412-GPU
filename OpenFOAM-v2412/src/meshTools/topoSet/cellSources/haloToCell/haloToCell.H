/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::haloToCell

Description
    A \c topoSetCellSource to select cells attached to the outside
    of this \c cellSet, and add into/remove from this \c cellSet.

    Operands:
    \table
      Operand  | Type    | Location
      input    | cellSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output   | cellSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      haloToCell;

        // Optional entries
        steps       1;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: haloToCell             | word |  yes  | -
      steps      | Number of grow/shrink steps to use  | label | no   | 1
    \endtable

    Options for the \c action entry:
    \verbatim
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

Note
    - \c action=new is not available for \c haloToCell.
    - When removing halo cells, this strips off
    any cells on the boundary of the cellSet.

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    haloToCell.C

\*---------------------------------------------------------------------------*/

#ifndef haloToCell_H
#define haloToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class haloToCell Declaration
\*---------------------------------------------------------------------------*/

class haloToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Number of grow/shrink steps
        label steps_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("haloToCell");


    // Constructors

        //- Construct from components
        explicit haloToCell(const polyMesh& mesh, const label nsteps=1);

        //- Construct from dictionary
        haloToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        haloToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~haloToCell() = default;


    // Member Functions

        //- The number of steps to grow/shrink
        label steps() const noexcept;

        //- Set number of steps to grow/shrink, return the old value
        label steps(const label nsteps) noexcept;


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

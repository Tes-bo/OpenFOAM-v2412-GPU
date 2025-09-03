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
    Foam::rotatedBoxToCell

Description
    A \c topoSetCellSource to select cells based on cell centres
    inside a given parallopiped (i.e. rotated/skewed box).

    Operands:
    \table
      Operand   | Type     | Location
      output    | cellSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      rotatedBoxToCell;
        origin      (<o1> <o2> <o3>);
        i           (<i1> <i2> <i3>);
        j           (<j1> <j2> <j3>);
        k           (<k1> <k2> <k3>);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: rotateTedBoxToCell     | word |  yes  | -
      origin     | Origin of the box                   | vector | yes | -
      i          | x1-axis                             | vector | yes | -
      j          | x2-axis                             | vector | yes | -
      k          | x3-axis                             | vector | yes | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Box defined as origin and i,j,k vectors.
    For example, box rotated 45 degrees around z-axis with
    sizes sqrt(0.2^2+0.2^2) (and extra large, 200 in z direction):
    \verbatim
      origin   ( 0.4 0.4 -100);
      i        ( 0.2 0.2  0);
      j        (-0.2 0.2  0);
      k        ( 0.0 0.0  200);
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    rotatedBoxToCell.C

\*---------------------------------------------------------------------------*/

#ifndef rotatedBoxToCell_H
#define rotatedBoxToCell_H

#include "topoSetCellSource.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class rotatedBoxToCell Declaration
\*---------------------------------------------------------------------------*/

class rotatedBoxToCell
:
    public topoSetCellSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Definition of the parallelopiped (i.e. rotated/skewed box)
        const point origin_;
        const vector i_;
        const vector j_;
        const vector k_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("rotatedBoxToCell");


    // Constructors

        //- Construct from components
        rotatedBoxToCell
        (
            const polyMesh& mesh,
            const vector& origin,
            const vector& i,
            const vector& j,
            const vector& k
        );

        //- Construct from dictionary
        rotatedBoxToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        rotatedBoxToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~rotatedBoxToCell() = default;


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

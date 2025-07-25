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
    Foam::clipPlaneToPoint

Description
    A \c topoSetPointSource to select all points
    above the specified plane

    Operands:
    \table
      Operand   | Type     | Location
      output    | pointSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      clipPlaneToPoint;

        point  (x y z);
        normal (x y z);
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of pointSet                    | word |  yes  | -
      type       | Type name: pointSet                 | word |  yes  | -
      action     | Action applied on points - see below| word |  yes  | -
      source     | Source name: clipPlaneToPoint        | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointSet from selected points
      add      | Add selected points into this pointSet
      subtract | Remove selected cells from this pointSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                        | Type   | Req'd  | Dflt
      point    | The reference point for the plane  | vector | yes    | -
      normal   | The unit normal of the plane       | vector | yes    | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    clipPlaneToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_clipPlaneToPoint_H
#define Foam_clipPlaneToPoint_H

#include "topoSetPointSource.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class clipPlaneToPoint Declaration
\*---------------------------------------------------------------------------*/

class clipPlaneToPoint
:
    public topoSetPointSource
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
    TypeName("clipPlaneToPoint");


    // Constructors

        //- No default construct
        clipPlaneToPoint() = delete;

        //- Construct from components
        clipPlaneToPoint
        (
            const polyMesh& mesh,
            const point& basePoint,
            const vector& normal
        );

        //- Construct from dictionary
        clipPlaneToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from stream - a single box.
        clipPlaneToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~clipPlaneToPoint() = default;


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

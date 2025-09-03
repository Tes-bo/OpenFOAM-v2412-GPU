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
    Foam::sphereToPoint

Description
    A \c topoSetPointSource to select points
    which are inside a given bounding sphere.

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
        source      sphereToPoint;
        origin      (0 0 0);
        radius      0.5;

        // Optional entries
        innerRadius 0;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                           | Type | Req'd | Dflt
      name       | Name of pointSet                      | word |  yes  | -
      type       | Type name: pointSet                   | word |  yes  | -
      action     | Action applied on points - see below  | word |  yes  | -
      source     | Source name: sphereToPoint            | word |  yes  | -
      origin     | The origin (centre) of the sphere     | vector | yes | -
      radius     | The (outer) radius of sphere          | scalar | yes | -
      innerRadius | The inner radius of sphere           | scalar | no  | 0
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
    sphereToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef sphereToPoint_H
#define sphereToPoint_H

#include "topoSetPointSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sphereToPoint Declaration
\*---------------------------------------------------------------------------*/

class sphereToPoint
:
    public topoSetPointSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Centre point of the sphere
        point origin_;

        //- The outer radius of the sphere
        scalar radius_;

        //- The inner radius of the sphere
        scalar innerRadius_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("sphereToPoint");


    // Constructors

        //- Construct from components
        sphereToPoint
        (
            const polyMesh& mesh,
            const point& origin,
            const scalar radius,
            const scalar innerRadius = 0
        );

        //- Construct from dictionary
        sphereToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        sphereToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~sphereToPoint() = default;


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

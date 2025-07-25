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
    Foam::labelToPoint

Description
    A \c topoSetPointSource to select points
    based on explicitly given point labels.

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
        source      labelToPoint;
        value
        (
            <pointLabel1>
            <pointLabel2>
            ...
        );
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Req'd | Dflt
      name       | Name of pointSet                     | word |  yes  | -
      type       | Type name: pointSet                  | word |  yes  | -
      action     | Action applied on points - see below | word |  yes  | -
      source     | Source name: labelToPoint            | word |  yes  | -
      value      | Input point labels                   | labelList | yes | -
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
    labelToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef labelToPoint_H
#define labelToPoint_H

#include "topoSetPointSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class labelToPoint Declaration
\*---------------------------------------------------------------------------*/

class labelToPoint
:
    public topoSetPointSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Point labels read from dictionary
        labelList labels_;


public:

    //- Runtime type information
    TypeName("labelToPoint");


    // Constructors

        //- Construct from components, copying labels
        labelToPoint(const polyMesh& mesh, const labelList& labels);

        //- Construct from components, moving labels
        labelToPoint(const polyMesh& mesh, labelList&& labels);

        //- Construct from dictionary
        labelToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        labelToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~labelToPoint() = default;


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

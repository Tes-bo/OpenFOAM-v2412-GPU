/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::setToPointZone

Description
    A \c topoSetSource to convert a \c pointSet
    to a \c pointZone (and associated \c pointSet).

    Operands:
    \table
      Operand   | Type      | Location
      input     | pointSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 1  | pointSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 2  | pointZone | $FOAM_CASE/constant/polyMesh/pointZones
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        pointZoneSet;
        action      <action>;

        // Mandatory entries
        source      setToCellZone;
        set         <pointSetName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                           | Type | Req'd | Dflt
      name       | Name of pointZone                     | word |  yes  | -
      type       | Type name: pointZoneSet               | word |  yes  | -
      action     | Action applied on points - see below  | word |  yes  | -
      source     | Source name: setToPointZone           | word |  yes  | -
      set        | Name of input pointSet                | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointZone from selected points of a pointSet
      add      | Add selected points of a pointSet into this pointZone
      subtract | Remove selected points of a pointSet from this pointZone
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    setToPointZone.C

\*---------------------------------------------------------------------------*/

#ifndef setToPointZone_H
#define setToPointZone_H

#include "topoSetPointZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setToPointZone Declaration
\*---------------------------------------------------------------------------*/

class setToPointZone
:
    public topoSetPointZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;


public:

    //- Runtime type information
    TypeName("setToPointZone");


    // Constructors

        //- Construct from components
        setToPointZone
        (
            const polyMesh& mesh,
            const word& setName
        );

        //- Construct from dictionary
        setToPointZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setToPointZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setToPointZone() = default;


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

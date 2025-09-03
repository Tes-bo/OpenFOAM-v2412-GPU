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
    Foam::setToCellZone

Description
    A \c topoSetSource to convert a \c cellSet
    to a \c cellZone (and associated \c cellSet).

    Operands:
    \table
      Operand   | Type     | Location
      input     | cellSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 1  | cellSet  | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output 2  | cellZone | $FOAM_CASE/constant/polyMesh/cellZones
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        cellZoneSet;
        action      <action>;

        // Mandatory entries
        source      setToCellZone;
        set         <cellSetName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellZone                    | word |  yes  | -
      type       | Type name: cellZoneSet              | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: setToCellZone          | word |  yes  | -
      set        | Name of input cellSet               | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellZone from selected cells of a cellSet
      add      | Add selected cells of a cellSet into this cellZone
      subtract | Remove selected cells of a cellSet from this cellZone
    \endverbatim

See also
    - Foam::topoSetSource

SourceFiles
    setToCellZone.C

\*---------------------------------------------------------------------------*/

#ifndef setToCellZone_H
#define setToCellZone_H

#include "topoSetCellZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setToCellZone Declaration
\*---------------------------------------------------------------------------*/

class setToCellZone
:
    public topoSetCellZoneSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of cellSet to use
        word setName_;


public:

    //- Runtime type information
    TypeName("setToCellZone");


    // Constructors

        //- Construct from components
        setToCellZone(const polyMesh& mesh, const word& setName);

        //- Construct from dictionary
        setToCellZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setToCellZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setToCellZone() = default;


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

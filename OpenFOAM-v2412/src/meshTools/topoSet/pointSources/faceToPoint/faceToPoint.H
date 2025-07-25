/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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
    Foam::faceToPoint

Description
    A \c topoSetPointSource to select all
    points based on usage in given \c faceSet(s).

    Operands:
    \table
      Operand   | Type       | Location
      input     | faceSet(s) | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | pointSet   | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      faceToPoint;
        option      <option>;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        sets
        (
            <faceSetName0>
            <faceSetName1>
            ...
        );

        // Option-2
        set     <faceSetName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Req'd | Dflt
      name       | Name of pointSet                     | word |  yes  | -
      type       | Type name: pointSet                  | word |  yes  | -
      action     | Action applied on points - see below | word |  yes  | -
      source     | Source name: faceToPoint             | word |  yes  | -
      option     | Selection type - see below           | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointSet from selected points
      add      | Add selected points into this pointSet
      subtract | Remove selected points from this pointSet
    \endverbatim

    Options for the \c option entry:
    \verbatim
      all       | Select all points of faces in the faceSet
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      sets     | Names of input faceSets    | wordList | cond'l | -
      zones    | Names of input faceZones   | wordList | cond'l | -
      set      | Name of input faceSet      | word     | cond'l | -
      zone     | Name of input faceZone     | word     | cond'l | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    faceToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef faceToPoint_H
#define faceToPoint_H

#include "topoSetPointSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceToPoint Declaration
\*---------------------------------------------------------------------------*/

class faceToPoint
:
    public topoSetPointSource
{

public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            ALL
        };


private:

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<faceAction> faceActionNames_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;

        //- Option
        faceAction option_;


    // Private Member Functions

        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("faceToPoint");


    // Constructors

        //- Construct from components
        faceToPoint
        (
            const polyMesh& mesh,
            const word& setName,
            const faceAction option
        );

        //- Construct from dictionary
        faceToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceToPoint() = default;


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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::pointToFace

Description
    A \c topoSetFaceSource to select faces with any
    point or any edge within a given \c pointSet(s).

    Operands
    \table
      Operand   | Type        | Location
      input     | pointSet(s) | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | faceSet     | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceSet;
        action      <action>;

        // Mandatory entries
        source      pointToFace;
        option      <option>;

        // Conditional mandatory entries
        // Select one of the below

        // Option-1
        sets
        (
            <pointSetName1>
            <pointSetName2>
            ...
        );

        // Option-2
        zones
        (
            <pointZoneName0>
            <pointZoneName1>
            ...
        );

        // Option-3
        set     <pointSetName>;

        // Option-4
        zone    <pointZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: pointToFace            | word |  yes  | -
      option     | Selection type - see below          | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

    Options for the \c option entry:
    \verbatim
      any     | Any point in the pointSet
      edge    | Two consecutive points (an edge) in the pointSet
      all     | All points in the pointSet
    \endverbatim

    Options for the conditional mandatory entries (in order of precedence):
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      sets     | Names of input pointSets   | wordList | cond'l | -
      zones    | Names of input pointZones  | wordList | cond'l | -
      set      | Name of input pointSet     | word     | cond'l | -
      zone     | Name of input pointZone    | word     | cond'l | -
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource

SourceFiles
    pointToFace.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointToFace_H
#define Foam_pointToFace_H

#include "topoSetFaceSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointToFace Declaration
\*---------------------------------------------------------------------------*/

class pointToFace
:
    public topoSetFaceSource
{

public:

        //- Enumeration defining the valid options
        enum pointAction
        {
            ANY,
            ALL,
            EDGE
        };


private:

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<pointAction> pointActionNames_;

        //- Names of sets or zones to use
        wordList names_;

        //- Is name a set or a zone
        const bool isZone_;

        //- Option
        pointAction option_;


    // Private Member Functions

        //- Depending on option add to or delete from set.
        template<class Selector>
        void combineImpl
        (
            topoSet& set,
            const bool add,
            const Selector& pointLabels
        ) const;

        //- Depending on face to cell option add to or delete from faceSet.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("pointToFace");


    // Constructors

        //- Construct from components
        pointToFace
        (
            const polyMesh& mesh,
            const word& setName,
            const pointAction option
        );

        //- Construct from dictionary
        pointToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        pointToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~pointToFace() = default;


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

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::planeToFaceZone

Description
    A \c topoSetSource to select faces based on the adjacent cell centres
    spanning a given plane. The plane is defined by a point and normal vector.

    Operands:
    \table
      Operand   | Type     | Location
      output 1  | faceSet  | \<constant\>/polyMesh/sets/\<set\>
      output 2  | faceZone | \<constant\>/polyMesh/faceZones
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        faceZoneSet;
        action      <action>;

        // Mandatory entries
        source      planeToFaceZone;
        point       (<px> <py> <pz>);
        normal      (<nx> <ny> <nz>);

        // Optional entries
        option      <option>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Reqd  | Dflt
      name       | Name of faceZone                    | word |  yes  | -
      type       | Type name: faceZoneSet              | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: planeToFaceZone        | word |  yes  | -
      point      | A point on the input plane          | vector | yes | -
      normal     | Normal vector to the input plane    | vector | yes | -
      option     | Plane selection method - see below  | word   | no  | all
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceZone from selected faces of a faceZone
      add      | Add selected faces of a faceZone into this faceZone
      subtract | Remove selected faces of a faceZone from this faceZone
    \endverbatim

    Options for the \c option entry:
    \verbatim
      all      | Select all faces that meet the criteria
      closest  | Select faces that belong to the closest contiguous plane
    \endverbatim

Note
  - Additionally, an include entry can be specified. When omitted or set to
    "all", then all faces that meet the criteria are included in the set. When
    set to "closest", just the faces that belong to the closest contiguous
    region to the plane point are included. This latter setting is useful when
    defining face zones through channels on which the flow rate is to be
    computed, as it keeps the set local to a single channel.

See also
    - Foam::topoSetSource

SourceFiles
    planeToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef planeToFaceZone_H
#define planeToFaceZone_H

#include "topoSetFaceZoneSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class planeToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class planeToFaceZone
:
    public topoSetFaceZoneSource
{
public:

        //- Enumeration defining the valid options
        enum faceAction
        {
            ALL,       //!< Select all faces that meet the criteria
            CLOSEST    //!< Select faces belong to the closest contiguous plane
        };


private:

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- Plane selection methods
        static const Enum<faceAction> faceActionNames_;

        //- Point on the input plane
        const vector point_;

        //- Normal to the input plane
        const vector normal_;

        //- Input plane selection method
        const faceAction option_;


    // Private Member Functions

        void combine(faceZoneSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("planeToFaceZone");


    // Constructors

        //- No default construct
        planeToFaceZone() = delete;

        //- Construct from components
        planeToFaceZone
        (
            const polyMesh& mesh,
            const point& basePoint,
            const vector& normal,
            const faceAction action = faceAction::ALL
        );

        //- Construct from dictionary
        planeToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        planeToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~planeToFaceZone() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

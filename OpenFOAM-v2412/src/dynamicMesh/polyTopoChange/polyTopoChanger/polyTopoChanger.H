/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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
    Foam::polyTopoChanger

Description
    List of mesh modifiers defining the mesh dynamics.

SourceFiles
    polyTopoChanger.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_polyTopoChanger_H
#define Foam_polyTopoChanger_H

#include "regIOobject.H"
#include "PtrList.H"
#include "polyMeshModifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class mapPolyMesh;
class polyBoundaryMesh;
class polyTopoChanger;

Ostream& operator<<(Ostream&, const polyTopoChanger&);


/*---------------------------------------------------------------------------*\
                      Class polyTopoChanger Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChanger
:
    public PtrList<polyMeshModifier>,
    public regIOobject
{
    // Private Member Functions

        //- Read if IOobject flags set, set modifiers. Return true if read.
        bool readIOcontents();

        //- No copy construct
        polyTopoChanger(const polyTopoChanger&) = delete;

        //- No copy assignment
        void operator=(const polyTopoChanger&) = delete;


protected:

    // Protected data

        //- Reference to mesh
        polyMesh& mesh_;

public:

    //- Runtime type information
    TypeName("polyTopoChanger");


    // Constructors

        //- Read construct given IOobject and a polyMesh
        polyTopoChanger(const IOobject& io, polyMesh& mesh);

        //- Read construct for given polyMesh and read-option
        polyTopoChanger(polyMesh& mesh, IOobjectOption::readOption rOpt);

        //- Read construct for given polyMesh.
        //  Uses read-option READ_IF_PRESENT
        explicit polyTopoChanger(polyMesh& mesh);


    //- Destructor
    virtual ~polyTopoChanger() = default;


    // Member functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of patch names
        wordList names() const;

        //- Is topology change required
        bool changeTopology() const;

        //- Return topology change request
        autoPtr<polyTopoChange> topoChangeRequest() const;

        //- Modify point motion
        void modifyMotionPoints(pointField&) const;

        autoPtr<mapPolyMesh> changeMesh
        (
            const bool inflate,
            const bool syncParallel = true,
            const bool orderCells = false,
            const bool orderPoints = false
        );

        //- Force recalculation of locally stored data on topological change
        void update(const mapPolyMesh& m);

        //- Add given set of topology modifiers to the topoChanger
        void addTopologyModifiers(const List<polyMeshModifier*>& tm);

        //- Find modifier given a name
        label findModifierID(const word& modName) const;


        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Member Operators

        bool operator!=(const polyTopoChanger&) const;
        bool operator==(const polyTopoChanger&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const polyTopoChanger&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

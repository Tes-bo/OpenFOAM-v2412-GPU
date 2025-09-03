/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::faceZone

Description
    A subset of mesh faces organised as a primitive patch.

    For quick check whether a face belongs to the zone use the lookup
    mechanism in faceZoneMesh, where all the zoned faces are registered
    with their zone number.

SourceFiles
    faceZone.C
    faceZoneNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faceZone_H
#define Foam_faceZone_H

#include "zone.H"
#include "faceZoneMeshFwd.H"
#include "primitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;
class faceZone;
Ostream& operator<<(Ostream& os, const faceZone& zn);


/*---------------------------------------------------------------------------*\
                          Class faceZone Declaration
\*---------------------------------------------------------------------------*/

class faceZone
:
    public zone
{
    // // Public Data Types
    //
    //     //- Side of the face zone
    //     enum side
    //     {
    //         FRONT = 1,  //!< The front (positive normal) side of the face
    //         BACK  = -1, //!< The back (negative normal) side of the face
    //     };


    // Private Data

        //- Flip map for all faces in the zone.
        //  True if the face needs flipping for the correct orientation.
        boolList flipMap_;

        //- Reference to zone list
        const faceZoneMesh& zoneMesh_;

        //- Demand-driven: Primitive patch of correctly flipped faces
        mutable std::unique_ptr<primitiveFacePatch> patchPtr_;

        //- Demand-driven: front cell layer (positive normal) side of faces
        mutable std::unique_ptr<labelList> frontCellsPtr_;

        //- Demand-driven: back cell layer (negative normal) side of faces
        mutable std::unique_ptr<labelList> backCellsPtr_;

        //- Demand-driven: Global edge addressing
        mutable std::unique_ptr<labelList> mePtr_;


    // Private Member Functions

        //- Set flip-map to uniform value
        void setFlipMap(const bool val);

        //- Build primitive patch
        void calcFaceZonePatch() const;

        //- Calculate front/back cell layers
        void calcCellLayers() const;

        //- Check addressing
        void checkAddressing() const;

        //- Clear out geometry and demand-driven data
        void clearGeom();

public:

    // Static Data Members

        //- The name associated with the zone-labels dictionary entry
        //- ("faceLabels")
        static const char * const labelsName;


    //- Runtime type information
    TypeName("faceZone");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceZone,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const label index,
                const faceZoneMesh& zm
            ),
            (name, dict, index, zm)
        );


    // Constructors

        //- No copy construct
        faceZone(const faceZone&) = delete;

        //- Construct an empty zone - name="", index=0
        explicit faceZone(const faceZoneMesh& zm);

        //- Construct an empty zone with specified name and index
        faceZone
        (
            const word& name,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components with uniform flip map value
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const bool flipMapValue,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components with uniform flip map value,
        //- transferring addressing.
        faceZone
        (
            const word& name,
            labelList&& addr,
            const bool flipMapValue,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components, copying addressing
        faceZone
        (
            const word& name,
            const labelUList& addr,
            const boolUList& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from components, transferring addressing
        faceZone
        (
            const word& name,
            labelList&& addr,
            boolList&& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct from dictionary
        faceZone
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct empty with original zone information (name, index, groups)
        //- and mesh reference. Optionally specify a new index.
        faceZone
        (
            const faceZone& originalZone,
            const Foam::zero,
            const faceZoneMesh& zm,
            const label newIndex = -1
        );

        //- Construct empty with original zone information (name, groups),
        //- resetting the index and zone mesh reference.
        faceZone
        (
            const faceZone& originalZone,
            const Foam::zero,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct with a new index and zone mesh information, the name
        //- of the original zone, resetting the face addressing
        //- and flip-map.
        faceZone
        (
            const faceZone& originalZone,
            const labelUList& addr,
            const boolUList& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct with a new index and zone mesh information, the name
        //- of the original zone, (move) resetting the face addressing
        //- and flip-map.
        faceZone
        (
            const faceZone& originalZone,
            labelList&& addr,
            boolList&& fm,
            const label index,
            const faceZoneMesh& zm
        );

        //- Construct and return a clone, resetting the zone mesh
        virtual autoPtr<faceZone> clone(const faceZoneMesh& zm) const
        {
            return autoPtr<faceZone>::New(*this, *this, flipMap(), index(), zm);
        }

        //- Construct and return a clone,
        //- resetting the face list and zone mesh
        virtual autoPtr<faceZone> clone
        (
            const labelUList& addr,
            const boolUList& fm,
            const label index,
            const faceZoneMesh& zm
        ) const
        {
            return autoPtr<faceZone>::New(*this, addr, fm, index, zm);
        }


    // Selectors

        //- Return a pointer to a new face zone
        //- created on freestore from dictionary
        static autoPtr<faceZone> New
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faceZoneMesh& zm
        );


    //- Destructor
    virtual ~faceZone() = default;


    // Member Functions

        //- Return reference to the zone mesh
        const faceZoneMesh& zoneMesh() const noexcept { return zoneMesh_; }

        //- The addressing (face IDs) used for the zone
        const labelList& addressing() const noexcept
        {
            return static_cast<const labelList&>(*this);
        }

        //- Return face flip map
        const boolList& flipMap() const noexcept { return flipMap_; }

        //- Helper function to re-direct to zone::localID(...)
        label whichFace(const label globalCellID) const;

        //- Return [demand-driven] reference to an equivalent primitive patch,
        //- with faces oriented according to flipMap()
        const primitiveFacePatch& patch() const;


    // Addressing into mesh

        //- The front cells layer.
        //- Cells on the positive normal side of faces.
        //  \warning may contain negative indices for boundary faces!
        const labelList& frontCells() const;

        //- The back cells layer.
        //- Cells on the negative normal side of faces.
        //  \warning may contain negative indices for boundary faces!
        const labelList& backCells() const;

        //- Return global edge index for local edges
        const labelList& meshEdges() const;

        //- Check zone definition. Return true if in error.
        virtual bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have faces synchronised.
        //  \return True if any errors.
        virtual bool checkParallelSync(const bool report = false) const;

        //- Correct patch after moving points
        virtual void movePoints(const pointField& pts);

        //- Update for changes in topology
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Write
        virtual void write(Ostream& os) const;

        //- Write dictionary
        virtual void writeDict(Ostream& os) const;


    // Assign addressing

        //- Clear addressing
        //- (remove lookup maps, patch/geometric information)
        virtual void clearAddressing();

        //- Clear primitive addressing
        virtual void clearPrimitives();

        //- Move reset addressing and flip map from another zone
        virtual void resetAddressing(faceZone&& zn);

        //- Copy reset addressing and flip map from another zone
        virtual void resetAddressing(const faceZone& zn);

        //- Reset addressing - use uniform flip map value
        //  Clears demand-driven data.
        virtual void resetAddressing
        (
            const labelUList& addr,
            const bool flipMapValue
        );

        //- Move reset addressing - use uniform flip map value
        //  Clears demand-driven data.
        virtual void resetAddressing
        (
            labelList&& addr,
            const bool flipMapValue
        );

        //- Copy reset addressing and flip map.
        //  Clears demand-driven data.
        virtual void resetAddressing
        (
            const labelUList& addr,
            const boolUList& flipMap
        );

        //- Assign addressing, clearing demand-driven data
        void operator=(const faceZone& zn);


    // I-O

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const faceZone& zn);


    // Housekeeping

        //- Deprecated(2023-09) same as patch()
        //  \deprecated(2023-09) same as patch()
        const primitiveFacePatch& operator()() const { return this->patch(); }

        //- Deprecated(2023-09) same as frontCells
        //  \deprecated(2023-09) same as frontCells #1832
        const labelList& masterCells() const { return frontCells(); }

        //- Deprecated(2023-09) same as backCells
        //  \deprecated(2023-09) same as backCells #1832
        const labelList& slaveCells() const { return backCells(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::wallDistAddressing

Description
    Variant of wallDist that uses meshWave and stores the addressing.

    Can be used to transport other wall properties.
    For example getting the nearest-wall normal into a volVectorField
    \verbatim
        // If not yet constructed construct from all wall patches
        const auto& wDist = wallDistAddressing::New(mesh);

        // Fill boundaryField with normal
        for (const label patchi : wDist.patchIDs())
        {
            auto& pnf = n.boundaryFieldRef()[patchi];
            pnf == pnf.patch().nf();
        }

        // Map data from nearest wall (using transformation if nearest is
        // reached through a coupled patch with transformation)
        wDist.map(n, mapDistribute::transform());
    \endverbatim


    Another use is e.g. to get the distance to cyclic patches. This enforces
    registration under the supplied name.

    \verbatim
        // If not yet constructed construct from selected patches
        const labelList& patchIDs = ..

        const auto& wDist = wallDistAddressing::New
        (
            "myPatches",    // registration name
            mesh,
            patchIDs
        );

        // Fill boundaryField with normal
        for (const label patchi : wDist.patchIDs())
        {
            auto& pnf = n.boundaryFieldRef()[patchi];
            pnf == pnf.patch().nf();
        }

        // Map data from nearest patch face (using transformation if nearest is
        // reached through a coupled patch with transformation)
        wDist.map(n, mapDistribute::transform());
    \endverbatim

Note
    By default (correctWalls = true) all
    cells point-connected to a wall explicitly search for the nearest
    location on the point-surrouding wall faces (so override the wave
    behaviour). This gets only done for processor-local wall faces.

See also
    Foam::wallDist

SourceFiles
    wallDistAddressing.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallDistAddressing_H
#define Foam_wallDistAddressing_H

#include "MeshObject.H"
#include "FaceCellWave.H"
#include "cellDistFuncs.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistribute;
class wallPointAddressing;

/*---------------------------------------------------------------------------*\
                     Class wallDistAddressing Declaration
\*---------------------------------------------------------------------------*/

class wallDistAddressing
:
    public MeshObject<fvMesh, UpdateableMeshObject, wallDistAddressing>,
    public cellDistFuncs
{
    // Private Typedefs

        typedef MeshObject
        <
            fvMesh,
            UpdateableMeshObject,
            wallDistAddressing
        > MeshObject_type;


    //- Use fvMesh& instead of cellDistFuncs polyMesh&
    using MeshObject_type::mesh_;


protected:

    // Private Data

        //- Set of patch IDs
        const labelList patchIDs_;

        //- Name for the patch set, e.g. "wall"
        const word patchTypeName_;

        //- Update wall distance every updateInterval_ steps
        const label updateInterval_;

        //- Do accurate distance calculation for near-wall cells.
        const bool correctWalls_;

        //- Flag to indicate whether the wall distance requires updating
        bool requireUpdate_;

        //- Distance-to-wall field
        mutable volScalarField y_;

        //- Number of wall faces
        autoPtr<globalIndex> globalWallsPtr_;

        //- Map to pull wall face info to cell or boundary face
        autoPtr<mapDistribute> mapPtr_;

        //- Indices of cells followed by boundary faces
        labelList untransformedItems_;

        //- Corresponding slot in mapPtr distribution result
        labelList untransformedSlots_;

        //- Start of patches. Start of untransformedPatchStarts_[0] is end
        //  of cell information
        labelList untransformedPatchStarts_;

        // Corresponding data for transformed items
        labelList transformedItems_;
        labelList transformedSlots_;
        labelList transformedPatchStarts_;


    // Protected Member Functions

        //- Extract FaceCellWave data
        label getValues
        (
            const FaceCellWave<wallPointAddressing>& wave,
            const List<wallPointAddressing>& allCellInfo,
            const List<wallPointAddressing>& allFaceInfo,
            volScalarField& y
        ) const;

        //- Store nearest-data to cell or boundary face
        void addItem
        (
            const label item,
            const labelPair& data,
            label& untransformi,
            label& transformi,
            labelPairList& transformedWallInfo
        );

        //- Extract nearest-patch distance data
        void correct(volScalarField& y);


        //- No copy construct
        wallDistAddressing(const wallDistAddressing&) = delete;

        //- No copy assignment
        void operator=(const wallDistAddressing&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("wallDistAddressing");


    // Constructors

        //- Construct from mesh and near-wall behaviour
        explicit wallDistAddressing
        (
            const fvMesh& mesh,
            const bool correctWalls = true,
            const label updateInterval = 1
        );

        //- Construct from patch type name (= registration name),
        //- mesh and patch IDs.
        wallDistAddressing
        (
            const word& patchTypeName,
            const fvMesh& mesh,
            const labelList& patchIDs,
            const bool correctWalls = true,
            const label updateInterval = 1
        );


    //- Destructor
    virtual ~wallDistAddressing();


    // Member Functions

        //- Return the patchIDs
        const labelUList& patchIDs() const noexcept
        {
            return patchIDs_;
        }

        //- Return reference to cached distance-to-wall field. Unvisited
        //  elements set to GREAT
        const volScalarField& y() const noexcept
        {
            return y_;
        }

        //- Update the y-field when the mesh moves
        virtual bool movePoints();

        //- Update the y-field when the mesh changes
        virtual void updateMesh(const mapPolyMesh&);


        // Get values from nearest patches

            //- Collect patchFields from patchIDs into straight list
            template<class Container, class Type>
            tmp<Field<Type>> collectPatchFields(const Container& bfld) const;

            //- Take collected/distributed patch field and fill volField
            template<class VolField>
            void extract
            (
                const Field<typename VolField::value_type>& patchFld,
                VolField& fld
            ) const;

            //- Map nearest-patch information. Take wall patch values
            //  and transports these to the internal field and other patches.
            //  Use mapDistribute::transformPosition if transporting absolute
            //  coordinates.
            template<class VolField, class TransformOp>
            const VolField& map
            (
                VolField& fld,
                const TransformOp& top = mapDistribute::transform()
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wallDistAddressingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020,2024 OpenCFD Ltd.
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
    Foam::cellDistFuncs

Description
    Collection of functions used in wall distance calculation.

SourceFiles
    cellDistFuncs.C
    cellDistFuncsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cellDistFuncs_H
#define cellDistFuncs_H

#include "HashSet.H"
#include "Map.H"
#include "wordRe.H"
#include "scalarField.H"
#include "point.H"
#include "primitivePatch.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class cellDistFuncs Declaration
\*---------------------------------------------------------------------------*/

class cellDistFuncs
{
    // Private Member Data

        //- Reference to mesh
        const polyMesh& mesh_;


    // Private Member Functions

        //- No copy construct
        cellDistFuncs(const cellDistFuncs&) = delete;

        //- No copy assignment
        void operator=(const cellDistFuncs&) = delete;

public:

    // Static Data Members

    ClassName("cellDistFuncs");


        //- Use combined-wall-patches wall distance v.s. v2406 per-patch
        //- distance. Default is true
        static bool useCombinedWallPatch;


    // Constructors

        //- Construct from mesh
        cellDistFuncs(const polyMesh& mesh);


    // Member Functions

        //- Access mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the set of patch IDs corresponding to the given names
        labelHashSet getPatchIDs(const UList<wordRe>& patchNames) const;

        //- Get patchIDs of/derived off certain type (e.g. 'processorPolyPatch')
        //  Uses isA, not isType
        template<class Type>
        labelHashSet getPatchIDs() const;

        //- Calculate smallest true distance (and patch face index)
        //  from pt to faces wallFaces.
        template<class PatchType>
        scalar smallestDist
        (
            const point& p,
            const PatchType& patch,
            const labelUList& wallFaces,
            label& patchFacei
        ) const;

        //- Get faces sharing point with face on patch
        template<class PatchType>
        void getPointNeighbours
        (
            const PatchType&,
            const label patchFacei,
            DynamicList<label>&
        ) const;

        //- Size of largest patch (out of supplied subset of patches)
        label maxPatchSize(const labelHashSet& patchIDs) const;

        //- Sum of patch sizes (out of supplied subset of patches).
        //  Used in sizing arrays.
        label sumPatchSize(const labelHashSet& patchIDs) const;

        //- Correct all cells connected to boundary (via face). Sets values in
        //  wallDistCorrected. Sets nearest wallface in nearestFace.
        void correctBoundaryFaceCells
        (
            const labelHashSet& patchIDs,
            scalarField& wallDistCorrected,
            Map<label>& nearestFace
        ) const;

        //- Correct all cells connected to wall (via point). Sets values in
        //  wallDistCorrected. Uses/sets nearest wallFace in nearestFace.
        void correctBoundaryPointCells
        (
            const labelHashSet& patchIDs,
            scalarField& wallDistCorrected,
            Map<label>& nearestFace
        ) const;

        //- Correct all cells connected to any of the patches in patchIDs. Sets
        //  - cell values in wallDistCorrected
        //  - (mesh) face that contains the nearest point
        void correctBoundaryCells
        (
            const labelList& patchIDs,
            const bool doPointCells,
            scalarField& wallDistCorrected,
            Map<label>& nearestFace
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cellDistFuncsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

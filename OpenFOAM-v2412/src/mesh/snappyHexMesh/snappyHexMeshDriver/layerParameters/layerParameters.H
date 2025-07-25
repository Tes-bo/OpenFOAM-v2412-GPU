/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::layerParameters

Description
    Simple container to keep together layer specific information.

SourceFiles
    layerParameters.C

\*---------------------------------------------------------------------------*/

#ifndef layerParameters_H
#define layerParameters_H

#include "dictionary.H"
#include "scalarField.H"
#include "labelList.H"
#include "Switch.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyBoundaryMesh;
class refinementSurfaces;

/*---------------------------------------------------------------------------*\
                           Class layerParameters Declaration
\*---------------------------------------------------------------------------*/

class layerParameters
{
public:

    // Public data types

        //- Enumeration defining the layer specification:
        //  - first and total thickness specified
        //  - first and expansion ratio specified
        //  - final and total thickness specified
        //  - final and expansion ratio specified
        //  - total thickness and expansion ratio specified
        enum thicknessModelType
        {
            FIRST_AND_TOTAL,
            FIRST_AND_EXPANSION,
            FINAL_AND_TOTAL,
            FINAL_AND_EXPANSION,
            TOTAL_AND_EXPANSION,
            FIRST_AND_RELATIVE_FINAL
        };


private:

    // Static data members

        //- Default angle for faces to be concave
        static const scalar defaultConcaveAngle;

        //- thicknessModelType names
        static const Enum<thicknessModelType> thicknessModelTypeNames_;


    // Private Data

        const dictionary dict_;

        //- In dry-run mode?
        const bool dryRun_;


        // Per patch (not region!) information

            //- How many layers to add.
            labelList numLayers_;

            //- Are sizes relative to local cell size
            boolList relativeSizes_;

            //- How thickness is specified.
            List<thicknessModelType> layerModels_;

            scalarField firstLayerThickness_;

            scalarField finalLayerThickness_;

            scalarField thickness_;

            scalarField expansionRatio_;

            //- Minimum total thickness
            scalarField minThickness_;


        const scalar featureAngle_;

        const scalar mergePatchFacesAngle_;

        const scalar concaveAngle_;

        const label nGrow_;

        const scalar maxFaceThicknessRatio_;

        const label nBufferCellsNoExtrude_;

        const label nLayerIter_;

        label nRelaxedIter_;

        //- Any additional reporting
        const bool additionalReporting_;

        const word meshShrinker_;

        const label nOuterIter_;


    // Private Member Functions

        //- Calculate expansion ratio from overall size v.s. thickness of
        //  first layer.
        static scalar layerExpansionRatio
        (
            const label n,
            const scalar totalOverFirst
        );

        //- Read parameters according to thickness model
        static void readLayerParameters
        (
            const bool verbose,
            const dictionary&,
            const thicknessModelType& spec,
            scalar& firstLayerThickness,
            scalar& finalLayerThickness,
            scalar& thickness,
            scalar& expansionRatio
        );

        void calculateLayerParameters
        (
            const thicknessModelType& spec,
            const label nLayers,
            scalar& firstLayerThickness,
            scalar& finalLayerThickness,
            scalar& thickness,
            scalar& expansionRatio
        );

        //- No copy construct
        layerParameters(const layerParameters&) = delete;

        //- No copy assignment
        void operator=(const layerParameters&) = delete;


public:

    // Constructors

        //- Construct from dictionary
        layerParameters
        (
            const dictionary& dict,
            const polyBoundaryMesh&,
            const bool dryRun = false
        );


    // Member Functions

        const dictionary& dict() const
        {
            return dict_;
        }


        // Per patch information

            //- How many layers to add.
            //  -1 : no specification. Assume 0 layers but allow sliding
            //       to make layers
            //   0 : specified to have 0 layers. No sliding allowed.
            //  >0 : number of layers
            const labelList& numLayers() const
            {
                return numLayers_;
            }

            //- Are size parameters relative to inner cell size or
            //  absolute distances.
            const boolList& relativeSizes() const
            {
                return relativeSizes_;
            }

            //- Specification of layer thickness
            const List<thicknessModelType>& layerModels() const
            {
                return layerModels_;
            }

            // Expansion factor for layer mesh
            const scalarField& expansionRatio() const
            {
                return expansionRatio_;
            }

            //- Wanted thickness of the layer furthest away
            //  from the wall (i.e. nearest the original mesh).
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& finalLayerThickness() const
            {
                return finalLayerThickness_;
            }

            //- Wanted thickness of the layer nearest to the wall.
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& firstLayerThickness() const
            {
                return firstLayerThickness_;
            }

            //- Wanted overall thickness of all layers.
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& thickness() const
            {
                return thickness_;
            }

            //- Minimum overall thickness of cell layer. If for any reason layer
            //  cannot be above minThickness do not add layer.
            //  If relativeSize() this number is relative to undistorted
            //  size of the cell outside layer.
            const scalarField& minThickness() const
            {
                return minThickness_;
            }


        // Control

            //- Number of overall layer addition iterations
            label nLayerIter() const
            {
                return nLayerIter_;
            }

            //- Outer loop to add layer by layer. Can be set to >= max layers
            //  in which case layers get added one at a time. This can help
            //  layer insertion since the newly added layers get included in
            //  the shrinking. Default is 1 -> add all layers in one go.
            label nOuterIter() const
            {
                return nOuterIter_;
            }

            //- Number of iterations after which relaxed motion rules
            //  are to be used.
            label nRelaxedIter() const
            {
                return nRelaxedIter_;
            }


        // Control

            scalar featureAngle() const
            {
                return featureAngle_;
            }

            scalar mergePatchFacesAngle() const
            {
                return mergePatchFacesAngle_;
            }

            scalar concaveAngle() const
            {
                return concaveAngle_;
            }

            //- If points get not extruded do nGrow layers of connected faces
            //  that are not grown. Is used to not do layers at all close to
            //  features.
            label nGrow() const
            {
                return nGrow_;
            }

            //- Stop layer growth on highly warped cells
            scalar maxFaceThicknessRatio() const
            {
                return maxFaceThicknessRatio_;
            }

            //- Create buffer region for new layer terminations
            label nBufferCellsNoExtrude() const
            {
                return nBufferCellsNoExtrude_;
            }

            //- Any additional reporting requested?
            bool additionalReporting() const
            {
                return additionalReporting_;
            }

            //- Type of mesh shrinker
            const word& meshShrinker() const
            {
                return meshShrinker_;
            }


        // Helper

            //- Determine overall thickness. Uses two of the four parameters
            //  according to the thicknessModel
            static scalar layerThickness
            (
                const thicknessModelType,
                const label nLayers,
                const scalar firstLayerThickness,
                const scalar finalLayerThickness,
                const scalar totalThickness,
                const scalar expansionRatio
            );

            //- Determine expansion ratio. Uses two of the four parameters
            //  according to the thicknessModel
            static scalar layerExpansionRatio
            (
                const thicknessModelType,
                const label nLayers,
                const scalar firstLayerThickness,
                const scalar finalLayerThickness,
                const scalar totalThickness,
                const scalar expansionRatio
            );

            //- Determine first layer (near-wall) thickness. Uses two of the
            //  four parameters according to the thicknessModel
            static scalar firstLayerThickness
            (
                const thicknessModelType,
                const label nLayers,
                const scalar firstLayerThickness,
                const scalar finalLayerThickness,
                const scalar totalThickness,
                const scalar expansionRatio
            );

            //- Determine ratio of final layer thickness to
            //  overall layer thickness
            static scalar finalLayerThicknessRatio
            (
                const label nLayers,
                const scalar expansionRatio
            );

            //- Determine overall thickness of a slice (usually 1 layer)
            static scalar layerThickness
            (
                const label nLayers,
                const scalar layerThickness,    // overall thickness
                const scalar expansionRatio,    // expansion ratio
                const label layerStart,         // starting layer (0=firstLayer)
                const label layerSize           // number of layers in slice
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

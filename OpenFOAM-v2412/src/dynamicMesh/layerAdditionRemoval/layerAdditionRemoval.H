/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::layerAdditionRemoval

Description
    Cell layer addition mesh modifier

SourceFiles
    layerAdditionRemoval.C
    addCellLayer.C
    removeCellLayer.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_layerAdditionRemoval_H
#define Foam_layerAdditionRemoval_H

#include "polyMeshModifier.H"
#include "primitivePatch.H"
#include "ZoneIDs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class layerAdditionRemoval Declaration
\*---------------------------------------------------------------------------*/

class layerAdditionRemoval
:
    public polyMeshModifier
{
    // Private Data

        //- Master face zone ID
        faceZoneID faceZoneID_;

        //- Min thickness of extrusion layer.  Triggers layer removal
        mutable scalar minLayerThickness_;

        //- Max thickness of extrusion layer.  Triggers layer addition
        mutable scalar maxLayerThickness_;

        //- Switch to calculate thickness as volume/area
        //  If false, thickness calculated from edges
        const bool thicknessFromVolume_;

        //- Layer thickness from previous step
        //  Used to decide whether to add or remove layers
        mutable scalar oldLayerThickness_;

        //- Point pairing
        mutable std::unique_ptr<labelList> pointsPairingPtr_;

        //- Face pairing
        mutable std::unique_ptr<labelList> facesPairingPtr_;

        //- Layer removal trigger time index
        mutable label triggerRemoval_;

        //- Layer addition trigger time index
        mutable label triggerAddition_;


    // Private Member Functions

        //- No copy construct
        layerAdditionRemoval(const layerAdditionRemoval&) = delete;

        //- No copy assignment
        void operator=(const layerAdditionRemoval&) = delete;

        //- Check validity of construction data
        void checkDefinition();


        // Topological Changes

            //- Check for valid layer
            bool validCollapse() const;

            //- Set layer pairing.  Return true if a valid layer exists
            bool setLayerPairing() const;

            //- Return points pairing in a layer (not automatic!)
            const labelList& pointsPairing() const;

            //- Return faces pairing in a layer (not automatic!)
            const labelList& facesPairing() const;

            //- Calculate the offset to the next layer
            tmp<vectorField> extrusionDir() const;

            //- Add a layer of cells
            void addCellLayer(polyTopoChange&) const;

            //- Remove a layer of cells
            void removeCellLayer(polyTopoChange&) const;

            //- Clear addressing
            void clearAddressing() const;


    // Static Data Members

        //- Thickness insertion fraction for the pre-motion
        static const scalar addDelta_;

        //- Thickness removal fraction for the cell collapse
        //  Note: the cell will be collapsed to this relative
        //  thickness before the layer is removed.
        static const scalar removeDelta_;

public:

    //- Runtime type information
    TypeName("layerAdditionRemoval");


    // Constructors

        //- Construct from components
        layerAdditionRemoval
        (
            const word& name,
            const label index,
            const polyTopoChanger& ptc,
            const word& zoneName,
            const scalar minThickness,
            const scalar maxThickness,
            const bool thicknessFromVolume = true
        );

        //- Construct from dictionary
        layerAdditionRemoval
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& ptc
        );


    //- Destructor
    virtual ~layerAdditionRemoval() = default;


    // Member Functions

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);


        // Edit

            //- Return min layer thickness which triggers removal
            scalar minLayerThickness() const
            {
                return minLayerThickness_;
            }

            //- Set min layer thickness which triggers removal
           void setMinLayerThickness(const scalar t) const;

            //- Return max layer thickness which triggers removal
            scalar maxLayerThickness() const
            {
                return maxLayerThickness_;
            }

            //- Set max layer thickness which triggers removal
           void setMaxLayerThickness(const scalar t) const;


        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

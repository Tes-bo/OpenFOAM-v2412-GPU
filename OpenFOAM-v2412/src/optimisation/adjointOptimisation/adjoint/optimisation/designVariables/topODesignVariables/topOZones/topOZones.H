/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    topOZones

Description
    Class holding cells zones used in topology optimisation

SourceFiles
    topOZones.C

\*---------------------------------------------------------------------------*/

#ifndef topOZones_H
#define topOZones_H

#include "fvMesh.H"
#include "betaMax.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class topOZones Declaration
\*---------------------------------------------------------------------------*/

class topOZones
{

protected:

    // Protected data

        //- Mesh reference
        const fvMesh& mesh_;

        //- TopO design variables dictionary
        const dictionary dict_;

        //- IDs of cellZones holding cells with constant alpha values
        labelList fixedPorousIDs_;

        //- The constant alpha values of fixedPorousIDs_
        scalarList fixedPorousValues_;

        //- IDs of cellZones holding cells with zero alpha values
        labelList fixedZeroPorousIDs_;

        //- IDs of cellZones with cells that can change their alpha value
        //- throughout the optimisation
        labelList adjointPorousIDs_;

        //- ID of the cellZone holding the cells next to inlets & outlets
        label IOcellsID_;

        //- The multiplier of the Brinkman penalisation term
        autoPtr<betaMax> betaMaxPtr_;


    // Protected Member Functions

        //- Get zone IDs corresponding to a wordList, read from a dict.
        //  Avoid going through ZoneMesh.indices() since this practically sorts
        //  the IDs from smallest to largest while we need to keep them in the
        //  same order as that prescribed in the wordList
        labelList getZoneIDs(const dictionary& dict, const word& zoneGroup);

        //- Add a cellZone containing the cells next to IO patches
        void addIOcellsZone();


private:

    // Private Member Functions

        //- No copy construct
        topOZones(const topOZones&);

        //- No copy assignment
        void operator=(const topOZones&);


public:

    //- Runtime type information
    TypeName("topOZones");


    // Constructors

        //- Construct from components
        topOZones
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~topOZones() = default;


    // Member Functions

        //- Const reference to mesh
        inline const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Cell zone IDs with fixed porosity values
        inline const labelList& fixedPorousZoneIDs() const
        {
            return fixedPorousIDs_;
        }

        //- Values of alpha for each fixed porous zone
        inline const scalarList& fixedPorousValues() const
        {
            return fixedPorousValues_;
        }

        //- Cell zone IDs with fixed zero porosity values
        inline const labelList& fixedZeroPorousZoneIDs() const
        {
            return fixedZeroPorousIDs_;
        }

        //- Cell zone IDs in which porosity is allowed to change
        inline const labelList& adjointPorousZoneIDs() const
        {
            return adjointPorousIDs_;
        }

        //- Cells next to IO boundaries
        inline const cellZone& IOCells() const
        {
            return mesh_.cellZones()[IOcellsID_];
        }

        //- ID of the cellZone holding the IOcells
        inline label IOzoneID() const
        {
            return IOcellsID_;
        }

        //- Get betaMax
        inline scalar getBetaMax() const
        {
            return betaMaxPtr_().value();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

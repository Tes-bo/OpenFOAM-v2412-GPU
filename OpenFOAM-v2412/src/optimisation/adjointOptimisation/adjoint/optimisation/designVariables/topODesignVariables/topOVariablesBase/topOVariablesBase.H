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
    Foam::topOVariablesBase

Description
    Base class for all design variables related to topology optimisation (topO).
    Provides the lookup functionality necessary for all fvOptions associated
    with topO and access to topOZones.

SourceFiles
    topOVariablesBase.C

\*---------------------------------------------------------------------------*/

#ifndef topOVariablesBase_H
#define topOVariablesBase_H

#include "localIOdictionary.H"
#include "topOZones.H"
#include "topOInterpolationFunction.H"
#include "FaceCellWave.H"
#include "wallPointData.H"
#include "cutFaceIso.H"
#include "cutCellIso.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class topOVariablesBase Declaration
\*---------------------------------------------------------------------------*/

class topOVariablesBase
:
    public localIOdictionary
{
protected:

    // Protected data

        //- Cell zones useful for defining the constant and changing parts
        //- of the domain in topO
        topOZones zones_;

        //- Folder name holding the zero level-set iso-surface
        fileName isoSurfFolder_;

        //- Map between iso-surface and mesh faces (internal and boundary)
        Map<label> meshFaceToChangedFace_;

        //- Per cutting face, the changed faces owned by it
        //labelListList changedFacesPerCuttingFace_;

        //- Iso-surface points
        pointField surfPoints_;

        //- Iso-surface faces
        faceList surfFaces_;


    // Protected Member Functions

        //- Based on the current pseudo-distance field y, identity the zero
        //- level-set iso-surface, write it to files for post-processing and
        //- set the seeds for the distance calculation from this iso-surface
        void setCutInterfaceFaces
        (
            const volScalarField& y,
            labelList& changedFaces,
            List<wallPointData<label>>& changedFacesInfo
        ) const;

        //- Construct facesFaces for a given boundary face
        DynamicList<label> faceFaces(const label facei) const;

        //- Check whether the cutFace intersects the boundary of the initial
        //- domain and add fluid part of the intersected face to the isoline
        bool addCutBoundaryFaceToIsoline
        (
            const label facei,
            const cutFaceIso& cutFace,
            DynamicList<vector>& isoSurfPts,
            DynamicList<face>& isoSurfFaces,
            DynamicList<label>& zoneIDs,
            List<DynamicList<label>>& cuttingFacesPerMeshFace
        ) const;

        bool isDuplicatePoint
        (
            const label pointID,
            const vector& pointi,
            const DynamicList<label>& cuttingFaces,
            const DynamicList<point>& isoSurfPts,
            const DynamicList<face>& isoSurfFaces,
            labelList& uniquePointIDs
        ) const;

        //- Add the cutting face to the zero level-set iso-surface
        bool addCuttingFaceToIsoline
        (
            const DynamicList<point>& facePoints,
            const label nSerialPatches,
            const DynamicList<label>& cellCutFaces,
            const List<DynamicList<label>>& cuttingFacesPerMeshFace,
            DynamicList<vector>& isoSurfPts,
            DynamicList<face>& isoSurfFaces,
            DynamicList<label>& zoneIDs
        ) const;

        //- Check whether the boundary faces of the initial domain belong
        //- to the fluid part and add them to the surface describing the
        //- fluid domain.
        //  Used only to close the output surface, does not set distance seeds
        void addBoundaryFacesToIsoline
        (
            const pointScalarField& pointY,
            const Map<label>& addedFaces,
            const scalar isoValue,
            DynamicList<vector>& isoSurfPts,
            DynamicList<face>& isoSurfFaces,
            DynamicList<label>& zoneIDs,
            label& nChangedFaces,
            labelList& changedFaces,
            List<wallPointData<label>>& changedFacesInfo,
            labelList& changedFaceToCutFace,
            List<DynamicList<label>>& cuttingFacesPerMeshFace
        );

        //- Write the surface describing the fluid domain to stl and vtp files
        void writeSurfaceFiles
        (
            const pointField& pts,
            const faceList& faces,
            const labelList& zoneIds,
            const label nSerialPatches
        ) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topOVariablesBase(const topOVariablesBase&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const topOVariablesBase&) = delete;


public:

    //- Runtime type information
    TypeName("topOVariablesBase");

    // Constructors

        //- Construct from dictionary
        topOVariablesBase
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~topOVariablesBase() = default;


    // Member Functions

        //- Get topOZones
        inline const topOZones& getTopOZones() const
        {
            return zones_;
        }

        //- Get betaMax value
        inline scalar getBetaMax() const
        {
            return zones_.getBetaMax();
        }

       //inline const labelListList& changedFacesPerCuttingFace() const
       //{
       //    return changedFacesPerCuttingFace_;
       //}

        inline const pointField& surfacePoints() const
        {
            return surfPoints_;
        }

        inline const faceList& surfaceFaces() const
        {
            return surfFaces_;
        }

        //- Get field used for physical interpolations
        virtual const volScalarField& beta() const = 0;

        //- Populate source terms for the flow equations
        virtual void sourceTerm
        (
            DimensionedField<scalar, volMesh>& field,
            const topOInterpolationFunction& interpolationFunc,
            const scalar betaMax,
            const word& interpolationFieldName = "beta"
        ) const;

        //- Post-processing sensitivities due to interpolations based on
        //- the indicator fields
        virtual void sourceTermSensitivities
        (
            scalarField& sens,
            const topOInterpolationFunction& interpolationFunc,
            const scalar betaMax,
            const word& designVariablesName,
            const word& interpolationFieldName = "beta"
        ) const;

        //- Write the fluid-solid interface to files.
        //  The interface is computed as an isoValue contour of the indicator
        //  field
        //  - 0 distance contour for levelSet or
        //  - 0.5 beta contour for porosity-based topO.
        //  For levelSet topO, the process of identifying the contour sets also
        //  the seeds for computing the distance field in the entire domain
        void writeFluidSolidInterface
        (
            const volScalarField& indicator,
            const scalar isoValue,
            labelList& changedFaces,
            List<wallPointData<label>>& changedFacesInfo
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

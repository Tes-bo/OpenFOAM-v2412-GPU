/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::sensitivitySurfacePoints

Description
    Calculation of adjoint-based sensitivities at wall points using the
    E-SI formulation

SourceFiles
    sensitivitySurfacePoints.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivitySurfacePoints_H
#define sensitivitySurfacePoints_H

#include "sensitivityShapeESI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class sensitivitySurfacePoints Declaration
\*---------------------------------------------------------------------------*/

class sensitivitySurfacePoints
:
    public sensitivityShapeESI
{
protected:

    // Protected data

        //- Write geometric info for use by external programs
        bool writeGeometricInfo_;

        //- Include surface area in sens computation
        bool includeSurfaceArea_;

        //- Is point belonging to a symmetry{Plane}
        boolList isSymmetryPoint_;

        //- Local point normal per symmetry point
        vectorField symmPointNormal_;;

        //- Extended patchIDs
        //  Sensitivities from patches adjacent to the sensitivityPatchIDs_
        //  should also be taken into consideration in order to compute the
        //  correct values at points in their interfaces
        labelHashSet extendedPatchIDs_;


    // Protected Member Functions

        //- Set suffix name for sensitivity fields
        labelHashSet populateExtendedIDs() const;

        //- Converts face sensitivities to point sensitivities and adds the
        //- ones directly computed in points (i.e. dSf/db and dnf/db).
        void finalisePointSensitivities();

        //- Construct globally correct point normals and point areas
        void constructGlobalPointNormalsAndAreas
        (
            vectorField& pointNormals,
            scalarField& pointMagSf
        );

        //- Set suffix name for sensitivity fields
        void setSuffixName();

        //- Allocate the proper size for the point-based sensitivities
        void computePointDerivativesSize();


private:

    // Private Member Functions

        //- No copy construct
        sensitivitySurfacePoints(const sensitivitySurfacePoints&) = delete;

        //- No copy assignment
        void operator=(const sensitivitySurfacePoints&) = delete;


public:

    //- Runtime type information
    TypeName("surfacePoints");


    // Constructors

        //- Construct from components
        sensitivitySurfacePoints
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver
        );


    //- Destructor
    virtual ~sensitivitySurfacePoints() = default;


    // Member Functions

        //- Read controls and update solver pointers if necessary
        void read();

        //- Read dict if changed
        virtual bool readDict(const dictionary& dict);

        //- Return set of patches on which to compute direct sensitivities
        virtual const labelHashSet& geometryVariationIntegrationPatches() const;

        //- Assemble sensitivities
        virtual void assembleSensitivities
        (
            autoPtr<designVariables>& designVars
        );

        //- Write sensitivity fields.
        //  If valid, copies boundaryFields to pointFields and writes them.
        virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

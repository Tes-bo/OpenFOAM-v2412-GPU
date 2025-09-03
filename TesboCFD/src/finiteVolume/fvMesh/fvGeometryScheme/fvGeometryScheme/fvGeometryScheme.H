/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::fvGeometryScheme

Description
    Abstract base class for geometry calculation schemes.

SourceFiles
    fvGeometryScheme.C

\*---------------------------------------------------------------------------*/

#ifndef fvGeometryScheme_H
#define fvGeometryScheme_H

#include "tmp.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class fvGeometryScheme Declaration
\*---------------------------------------------------------------------------*/

class fvGeometryScheme
:
    public refCount
{
    // Private Member Functions

        //- No copy construct
        fvGeometryScheme(const fvGeometryScheme&) = delete;

        //- No copy assignment
        void operator=(const fvGeometryScheme&) = delete;


protected:

    // Protected Data

        //- Hold reference to mesh
        const fvMesh& mesh_;


    // Protected Member Functions

        //- Set the mesh motion flux
        bool setMeshPhi() const;


public:

    //- Runtime type information
    TypeName("fvGeometryScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            fvGeometryScheme,
            dict,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from mesh
        fvGeometryScheme(const fvMesh& mesh, const dictionary& dict)
        :
            mesh_(mesh)
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<fvGeometryScheme> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& defaultScheme
        );


    //- Destructor
    virtual ~fvGeometryScheme() = default;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Update basic geometric properties from provided points
        virtual void movePoints();

        //- Update mesh for topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Return linear difference weighting factors
        virtual tmp<surfaceScalarField> weights() const = 0;

        //- Return cell-centre difference coefficients
        virtual tmp<surfaceScalarField> deltaCoeffs() const = 0;

        //- Return non-orthogonal cell-centre difference coefficients
        virtual tmp<surfaceScalarField> nonOrthDeltaCoeffs() const = 0;

        //- Return non-orthogonality correction vectors
        virtual tmp<surfaceVectorField> nonOrthCorrectionVectors() const = 0;

        ////- Selector for wall distance method. WIP. Ideally return wall
        ////  distance or meshObject?
        //virtual autoPtr<patchDistMethod> newPatchDistMethod
        //(
        //    const dictionary& dict,
        //    const labelHashSet& patchIDs,
        //    const word& defaultPatchDistMethod
        //) const = 0;

        //- Calculate geometry quantities using mesh topology and provided
        //- points. If oldPoints provided only does local update. Returns
        //- true if anything changed, false otherwise
        virtual bool updateGeom
        (
            const pointField& points,
            const refPtr<pointField>& oldPoints,    // optional old points
            pointField& faceCentres,
            vectorField& faceAreas,
            pointField& cellCentres,
            scalarField& cellVolumes
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

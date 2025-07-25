/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::deltaBoundary

Description
    Differentiation of the mesh data structure

SourceFiles
    deltaBoundary.C

\*---------------------------------------------------------------------------*/

#ifndef deltaBoundary_H
#define deltaBoundary_H

#include "fieldTypes.H"
#include "vectorField.H"
#include "tensorField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class deltaBoundary Declaration
\*---------------------------------------------------------------------------*/

class deltaBoundary
{
protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;


private:

    // Private Member Functions

        //- No copy construct
        deltaBoundary(const deltaBoundary&) = delete;

        //- No copy assignment
        void operator=(const deltaBoundary&) = delete;

        //- Compute tensor-vector products
        tensor tensorCrossVector(const tensor& T, const vector& v);


public:

    // Constructors

        //- Construct from components
        deltaBoundary(const fvMesh& mesh);


    //- Destructor
    ~deltaBoundary() = default;


    // Member Functions

        //- Given a face and the points to be moved in the normal direction,
        //- find faceArea, faceCentre and unitVector changes
        vectorField makeFaceCentresAndAreas_d
        (
            const pointField& p,
            const pointField& p_d
        );

        //- Given a face and the points to be moved in an arbitrary direction,
        //- find faceArea, faceCentre and unitVector changes
        tensorField makeFaceCentresAndAreas_d
        (
            const pointField& p,
            const tensorField& p_d
        );

        //- Compute cell center variation wrt
        //- given face movement or derivative.
        //  pT (perturbation type) should be a vector in case
        //  of known face area and ctr movements
        //  or a tensor for gradients
        template<class pT>
        pT makeCellCentres_d
        (
            const vectorField& fAreas,
            const vectorField& fCtrs,
            const Field<pT>& fAreas_d,
            const Field<pT>& fCtrs_d
        );

        //- Compute the change of the cell centers of the pointCells
        //- of pointI, for a unitary movement of pointI in all three directions
        tmp<tensorField> cellCenters_d(const label pointI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "deltaBoundaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

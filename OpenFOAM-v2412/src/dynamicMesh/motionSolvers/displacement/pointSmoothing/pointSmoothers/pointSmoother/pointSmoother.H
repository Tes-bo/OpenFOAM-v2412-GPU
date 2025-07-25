/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::pointSmoother

Description
    Abstract base class for point smoothing methods. Handles parallel
    communication via reset and average functions.

SourceFiles
    pointSmoother.C
    pointSmootherTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef pointSmoother_H
#define pointSmoother_H

#include "polyMeshGeometry.H"
#include "runTimeSelectionTables.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointSmoother Declaration
\*---------------------------------------------------------------------------*/

class pointSmoother
{
private:

    // Private data

        //- Reference to the polyMesh
        const polyMesh& mesh_;

        //- Set of the processor patch indices
        labelHashSet processorPatchIDs_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointSmoother(const pointSmoother&);

        //- Disallow default bitwise assignment
        void operator=(const pointSmoother&);


protected:

    // Protected Member Functions

        //- Test if the given face is internal or on a processor boundary
        bool isInternalOrProcessorFace(const label faceI) const;

        //- Get a boolean list of the points to be moved
        autoPtr<PackedBoolList> pointsToMove
        (
            const labelList& facesToMove,
            const bool moveInternalFaces
        ) const;

        //- Reset the relevant weights and displacements to zero
        template <class weightType>
        void reset
        (
            const labelList& facesToMove,
            Field<weightType>& weights,
            vectorField& pointDisplacement,
            const bool resetInternalFaces = true
        ) const;

        //- Average the displacements using the weights provided
        template <class weightType>
        void average
        (
            const labelList& facesToMove,
            Field<weightType>& weights,
            vectorField& pointDisplacement
        ) const;


public:

    //- Runtime type information
    TypeName("pointSmoother");


    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        pointSmoother,
        dictionary,
        (const polyMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );


    // Constructors

        //- Construct from a dictionary and a point displacement field
        pointSmoother(const polyMesh& mesh, const dictionary&);


    // Selector

        //- Construct given type
        static autoPtr<pointSmoother> New
        (
            const word& pointSmootherType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct with type looked up from dictionary
        static autoPtr<pointSmoother> New
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~pointSmoother();


    // Member Functions

        // Access the mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Update the point displacements and apply constraints
        void update
        (
            const labelList& facesToMove,
            const pointField& oldPoints,
            const pointField& currentPoints,
            const pointField& faceCentres,
            const vectorField& faceAreas,
            const pointField& cellCentres,
            const scalarField& cellVolumes,
            pointVectorField& pointDisplacement,
            const bool correctBCs = true
        ) const;

        //- Update the point displacements and apply constraints
        void update
        (
            const labelList& facesToMove,
            const pointField& oldPoints,
            const pointField& currentPoints,
            const polyMeshGeometry& meshGeometry,
            pointVectorField& pointDisplacement,
            const bool correctBCs = true
        ) const
        {
            update
            (
                facesToMove,
                oldPoints,
                currentPoints,
                meshGeometry.faceCentres(),
                meshGeometry.faceAreas(),
                meshGeometry.cellCentres(),
                meshGeometry.cellVolumes(),
                pointDisplacement,
                correctBCs
            );
        }

        //- Calculate the point displacement
        virtual void calculate
        (
            const labelList& facesToMove,
            const pointField& oldPoints,
            const pointField& currentPoints,
            const pointField& faceCentres,
            const vectorField& faceAreas,
            const pointField& cellCentres,
            const scalarField& cellVolumes,
            vectorField& pointDisplacement
        ) const = 0;

        //- Update the point displacements
        virtual void calculate
        (
            const labelList& facesToMove,
            const pointField& oldPoints,
            const pointField& currentPoints,
            const polyMeshGeometry& meshGeometry,
            vectorField& pointDisplacement
        ) const
        {
            calculate
            (
                facesToMove,
                oldPoints,
                currentPoints,
                meshGeometry.faceCentres(),
                meshGeometry.faceAreas(),
                meshGeometry.cellCentres(),
                meshGeometry.cellVolumes(),
                pointDisplacement
            );
        }

        //- Check element quality: 1 = best, 0 = invalid. (also negative?)
        //- Topology from mesh, point locations supplied.
        //- Move to motionSolver level?
        virtual tmp<scalarField> faceQuality
        (
            const pointField& points,
            const pointField& faceCentres,
            const vectorField& faceAreas,
            const pointField& cellCentres,
            const scalarField& cellVolumes
        ) const;

        //- Check element quality: 1 = best, 0 = invalid.
        //- Topology from mesh, point locations supplied.
        //- Move to motionSolver level?
        virtual tmp<scalarField> cellQuality
        (
            const pointField& points,
            const pointField& faceCentres,
            const vectorField& faceAreas,
            const pointField& cellCentres,
            const scalarField& cellVolumes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pointSmootherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

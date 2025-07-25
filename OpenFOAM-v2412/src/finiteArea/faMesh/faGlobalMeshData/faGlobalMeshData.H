/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::faGlobalMeshData

Description
    Various mesh related information for a parallel run

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef Foam_faGlobalMeshData_H
#define Foam_faGlobalMeshData_H

#include "processorTopology.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class Declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                        Class faGlobalMeshData Declaration
\*---------------------------------------------------------------------------*/

class faGlobalMeshData
{
    // Private Data

        //- Reference to mesh
        const faMesh& mesh_;

        //- The processor/processor topology
        processorTopology processorTopology_;


    // Globally shared point addressing

        //- Total number of global points
        label nGlobalPoints_;

        //- Indices of local points that are globally shared
        labelList sharedPointLabels_;

        //- Indices of globally shared points in the master list
        //  This list contains all the shared points in the mesh
        labelList sharedPointAddr_;


public:

    //- Runtime type information
    ClassName("faGlobalMeshData");


    // Generated Methods

        //- No copy construct
        faGlobalMeshData(const faGlobalMeshData&) = delete;

        //- No copy assignment
        void operator=(const faGlobalMeshData&) = delete;


    // Constructors

        //- Construct from mesh
        explicit faGlobalMeshData(const faMesh& mesh);


    //- Destructor
    ~faGlobalMeshData();


    // Member Functions

    // Access

        //- Return mesh reference
        const faMesh& mesh() const noexcept;


    // Processor-Topology

        //- The processor to processor topology.
        const processorTopology& topology() const noexcept
        {
            return processorTopology_;
        }

        //- Order in which the patches should be initialised/evaluated
        //- corresponding to the schedule
        const lduSchedule& patchSchedule() const noexcept
        {
            return processorTopology_.patchSchedule();
        }


    // Globally shared point addressing

        //- Return number of globally shared points
        label nGlobalPoints() const noexcept
        {
            return nGlobalPoints_;
        }

        //- Return indices of local points that are globally shared
        const labelList& sharedPointLabels() const noexcept
        {
            return sharedPointLabels_;
        }

        //- Return addressing into the complete globally shared points list
        const labelList& sharedPointAddr() const noexcept
        {
            return sharedPointAddr_;
        }

        //- Change global mesh data given a topological change.
        void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

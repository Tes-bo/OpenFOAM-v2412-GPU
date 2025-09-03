/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::parLagrangianDistributor

Description
    Lagrangian field redistributor.

    Runs in parallel. Redistributes from fromMesh to toMesh.

SourceFiles
    parLagrangianDistributor.C
    parLagrangianDistributorFields.C
    parLagrangianDistributorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_parLagrangianDistributor_H
#define Foam_parLagrangianDistributor_H

#include "PtrList.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistributePolyMesh;
class mapDistributeBase;
class IOobjectList;
class passivePositionParticleCloud;

/*---------------------------------------------------------------------------*\
                    Class parLagrangianDistributor Declaration
\*---------------------------------------------------------------------------*/

class parLagrangianDistributor
{
    // Private Data

        //- Source mesh reference
        const fvMesh& srcMesh_;

        //- Destination mesh reference
        const fvMesh& tgtMesh_;

        //- Distribution map reference
        const mapDistributePolyMesh& distMap_;

        //- For every src cell the target processor
        labelList destinationProcID_;

        //- For every src cell the target cell
        labelList destinationCell_;


public:

    //- Output verbosity when writing
    static int verbose_;


    // Generated Methods

        //- No copy construct
        parLagrangianDistributor(const parLagrangianDistributor&) = delete;

        //- No copy assignment
        void operator=(const parLagrangianDistributor&) = delete;


    // Constructors

        //- Construct from components
        parLagrangianDistributor
        (
            const fvMesh& srcMesh,
            const fvMesh& tgtMesh,
            const label nOldCells,
            const mapDistributePolyMesh& distMap
        );


    // Static Functions

        //- Find all clouds (on all processors) and for each cloud all
        //- the objects. Result will be synchronised on all processors
        static void findClouds
        (
            const fvMesh&,

            //! All cloud names on any processor
            wordList& cloudNames,

            //! Per cloud name, whether my processor has it
            boolList& haveClouds,

            //! Per cloud nmae, the field names
            List<wordList>& objectNames
        );

        //- Pick up any fields of a given type
        template<class Type>
        static wordList filterObjects
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Read and store all fields of a cloud
        template<class Container>
        static label readFields
        (
            const passivePositionParticleCloud& cloud,
            const bool haveCloud,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Read and store all fields for known cloud field types
        static label readAllFields
        (
            const passivePositionParticleCloud& cloud,
            const bool haveCloud,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );


    // Member Functions

        //- Redistribute and write lagrangian positions
        autoPtr<mapDistributeBase> distributeLagrangianPositions
        (
            passivePositionParticleCloud& cloud
        ) const;

        //- Read, redistribute and write lagrangian positions
        autoPtr<mapDistributeBase> distributeLagrangianPositions
        (
            const word& cloudName
        ) const;

        //- Redistribute all fields for known cloud field types
        label distributeAllFields
        (
            const mapDistributeBase& lagrangianMap,
            const word& cloudName,
            const bool haveCloud,
            const IOobjectList& cloudObjs,
            const wordRes& selectedFields
        ) const;

        //- Redistribute and write all stored lagrangian fields
        label distributeAllStoredFields
        (
            const mapDistributeBase& lagrangianMap,
            passivePositionParticleCloud& cloud
        ) const;


        //- Read, redistribute and write all/selected lagrangian fields
        template<class Type>
        label distributeFields
        (
            const mapDistributeBase& map,
            const word& cloudName,
            const bool haveCloud,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Read, redistribute and write all/selected lagrangian fieldFields
        template<class Type>
        label distributeFieldFields
        (
            const mapDistributeBase& map,
            const word& cloudName,
            const bool haveCloud,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Redistribute and write stored lagrangian fields.
        //  Note: does no reading so no need to check for existence
        //  of lagrangian files
        template<class Container>
        label distributeStoredFields
        (
            const mapDistributeBase& map,
            passivePositionParticleCloud& cloud
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "parLagrangianDistributorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

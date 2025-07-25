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

Namespace
    Foam::DMDModels

Description
    A namespace for various dynamic mode
    decomposition (DMD) model implementations.

Class
    Foam::DMDModel

Description
    Abstract base class for DMD models to handle DMD
    characteristics for the \c DMD function object.

SourceFiles
    DMDModel.C
    DMDModelNew.C
    DMDModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef DMDModel_H
#define DMDModel_H

#include "fvMesh.H"
#include "dictionary.H"
#include "HashSet.H"
#include "runTimeSelectionTables.H"
#include "writeFile.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class DMDModel Declaration
\*---------------------------------------------------------------------------*/

class DMDModel
:
    public functionObjects::writeFile
{
    typedef SquareMatrix<scalar> SMatrix;
    typedef RectangularMatrix<scalar> RMatrix;


protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Name of operand function object
        const word name_;


    // Protected Member Functions

        // Evaluation

            //- Compute and write mode dynamics
            virtual bool dynamics() = 0;

            //- Compute and write modes
            virtual bool modes() = 0;


public:

    //- Runtime type information
    TypeName("DMDModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            DMDModel,
            dictionary,
            (
                const fvMesh& mesh,
                const word& name,
                const dictionary& dict
            ),
            (mesh, name, dict)
        );


    // Selectors

        //- Return a reference to the selected DMD model
        static autoPtr<DMDModel> New
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        DMDModel
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& dict
        );

        //- No copy construct
        DMDModel(const DMDModel&) = delete;

        //- No copy assignment
        void operator=(const DMDModel&) = delete;


    //- Destructor
    virtual ~DMDModel() = default;


    // Member Functions

        // Evaluation

            //- Initialise model data members with a given snapshot
            virtual bool initialise(const RMatrix& snapshot) = 0;

            //- Update model data members with a given snapshot
            virtual bool update(const RMatrix& snapshot) = 0;

            //- Compute and write modes and
            //- mode dynamics of model data members
            virtual bool fit() = 0;

            //- Compute and write a reconstruction of flow field
            //- based on given modes and mode dynamics (currently no-op)
            virtual void reconstruct(const wordList modes)
            {
                NotImplemented;
            }


        // I-O

            //- Read model settings
            virtual bool read(const dictionary& dict) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

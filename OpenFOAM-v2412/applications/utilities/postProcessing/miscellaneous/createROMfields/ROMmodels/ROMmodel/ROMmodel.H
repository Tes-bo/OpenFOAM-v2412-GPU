/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::ROMmodels

Description
    A namespace for various implementations of
    reduced-order (ROM) field creation models.

Class
    Foam::ROMmodel

Description
    Abstract base class for reduced-order models
    to handle specific model characteristics.

SourceFiles
    ROMmodel.C
    ROMmodelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_ROMmodel_H
#define Foam_ROMmodel_H

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Time;
class fvMesh;
class dictionary;
class instant;
typedef class List<instant> instantList;

/*---------------------------------------------------------------------------*\
                          Class ROMmodel Declaration
\*---------------------------------------------------------------------------*/

class ROMmodel
{
protected:

    // Protected Data

        //- Reference to the Time
        //  Need non-const access to use setTime
        Time& runTime_;

        //- Reference to the fvMesh
        //  Need non-const access to use readFieldsHandler
        fvMesh& mesh_;

        //- Const reference to the dictionary
        const dictionary& dict_;

        //- Const reference to field times
        const instantList& times_;


public:

    //- Runtime type information
    TypeName("ROMmodel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ROMmodel,
            dictionary,
            (
                Time& runTime,
                fvMesh& mesh,
                const dictionary& dict,
                const instantList& times
            ),
            (runTime, mesh, dict, times)
        );


    // Selectors

        //- Return a reference to the selected ROMmodel
        static autoPtr<ROMmodel> New
        (
            Time& runTime,
            fvMesh& mesh,
            const dictionary& dict,
            const instantList& times
        );


    // Constructors

        //- Construct from components
        ROMmodel
        (
            Time& runTime,
            fvMesh& mesh,
            const dictionary& dict,
            const instantList& times
        );

        //- No copy construct
        ROMmodel(const ROMmodel&) = delete;

        //- No copy assignment
        void operator=(const ROMmodel&) = delete;


    //- Destructor
    virtual ~ROMmodel() = default;


    // Member Functions

        //- Read model settings
        virtual bool read(const dictionary& dict) = 0;

        //- Create and write fields
        virtual bool createAndWrite() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

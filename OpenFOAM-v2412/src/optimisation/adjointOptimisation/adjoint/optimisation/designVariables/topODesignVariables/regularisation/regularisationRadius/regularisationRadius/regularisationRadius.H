/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 PCOpt/NTUA
    Copyright (C) 2021 FOSS GP
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
    Foam::regularisationRadius

Description
    Base class for selecting the regulatisation radius


SourceFiles
    regularisationRadius.C

\*---------------------------------------------------------------------------*/

#ifndef regularisationRadius_H
#define regularisationRadius_H

#include "fvMatrices.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class regularisationRadius Declaration
\*---------------------------------------------------------------------------*/

class regularisationRadius
{

private:

    // Private Member Functions

        //- No copy construct
        regularisationRadius(const regularisationRadius&) = delete;

        //- No copy assignment
        void operator=(const regularisationRadius&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;

        const dictionary dict_;


public:

    //- Runtime type information
    TypeName("regularisationRadius");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            regularisationRadius,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                bool adjustWallThickness
            ),
            (mesh, dict, adjustWallThickness)
        );


    // Constructors

        //- Construct from components
        regularisationRadius
        (
            const fvMesh& mesh,
            const dictionary& dict,
            bool adjustWallThickness
        );


    // Selectors

        //- Construct and return the selected regularisationRadius
        //  Default to isotropic
        static autoPtr<regularisationRadius> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            bool adjustWallThickness
        );


    //- Destructor
    virtual ~regularisationRadius() = default;


    // Member Functions

        //- Term including the regulatisation of the field
        virtual void addRegularisationTerm
        (
            fvScalarMatrix& matrix,
            bool isTopoField
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

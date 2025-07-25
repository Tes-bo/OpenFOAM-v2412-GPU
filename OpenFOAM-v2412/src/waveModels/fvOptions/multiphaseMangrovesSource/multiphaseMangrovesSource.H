/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::fv::multiphaseMangrovesSource

Group
    grpFvOptionsSources

Description

Usage

SourceFiles
    multiphaseMangrovesSource.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseMangrovesSource_H
#define multiphaseMangrovesSource_H

#include "fvOption.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class MangrovesModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class multiphaseMangrovesSource Declaration
\*---------------------------------------------------------------------------*/

class multiphaseMangrovesSource
:
    public fv::option
{
    // Private Member Functions

        //- No copy construct
        multiphaseMangrovesSource(const multiphaseMangrovesSource&) = delete;

        //- No copy assignment
        void operator=(const multiphaseMangrovesSource&) = delete;


protected:

    // Protected data

        // Coefficients per cell zone

            //- Width of the vegetation element
            scalarList aZone_;

            //- Number of plants per unit of area
            scalarList NZone_;

            //- Inertia coefficient
            scalarList CmZone_;

            //- Drag coefficient
            scalarList CdZone_;

            //- Porosity coefficient
            //scalarList MangrovesZone_;

            //- Zone indices
            labelListList zoneIDs_;


    // Protected Member Functions

        //- Return the drag force coefficient
        tmp<volScalarField> dragCoeff(const volVectorField& U) const;

        //- Return the inertia force coefficient
        tmp<volScalarField> inertiaCoeff() const;


public:

    //- Runtime type information
    TypeName("multiphaseMangrovesSource");


    // Constructors

        //- Construct from components
        multiphaseMangrovesSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~multiphaseMangrovesSource() = default;


    // Member Functions

        // Add explicit and implicit contributions

            //- Add implicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add implicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

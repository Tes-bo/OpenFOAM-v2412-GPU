/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::regionModels::thermalShellModels::thermalShellModel

Description
    Intermediate class for thermal-shell finite-area models.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        T                     <word>;

        // Optional entries
        thermalShellModel     <word>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                       | Type  | Reqd | Deflt
      T          | Name of operand temperature field | word  | yes  | -
      thermalShellModel | Name of thermal-shell model <!--
                        -->                      | word | choice | -
    \endtable

    The inherited entries are elaborated in:
      - \link regionFaModel.H \endlink

SourceFiles
    thermalShellModel.C
    thermalShellModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_regionModels_thermalShellModel_H
#define Foam_regionModels_thermalShellModel_H

#include "areaFieldsFwd.H"
#include "volFieldsFwd.H"
#include "regionFaModel.H"
#include "faOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                      Class thermalShellModel Declaration
\*---------------------------------------------------------------------------*/

class thermalShellModel
:
    public regionFaModel
{
protected:

    // Protected Data

        //- Name of the temperature field
        word TName_;

        //- Primary region temperature
        const volScalarField& Tp_;

        //- Shell temperature
        areaScalarField T_;

        //- Reference to faOptions
        Foam::fa::options& faOptions_;


public:

    //- Runtime type information
    TypeName("thermalShellModel");


    // Declare runtime constructor selection tables

         declareRunTimeSelectionTable
         (
             autoPtr,
             thermalShellModel,
             dictionary,
             (
                const word& modelType,
                const fvMesh& mesh,
                const dictionary& dict
             ),
             (modelType, mesh, dict)
         );


    // Constructors

        //- Construct from type name and mesh and dict
        thermalShellModel
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        thermalShellModel(const thermalShellModel&) = delete;

        //- No copy assignment
        void operator=(const thermalShellModel&) = delete;


    // Selectors

        //- Return a reference to the selected model using dictionary
        static autoPtr<thermalShellModel> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~thermalShellModel() = default;


    // Member Functions

        // Access

            //- Return primary region temperature
            const volScalarField& Tp() const noexcept
            {
                return Tp_;
            }

            //- Return shell temperature
            const areaScalarField& T() const noexcept
            {
                return T_;
            }

            //- Return faOptions
            Foam::fa::options& faOptions() noexcept
            {
                return faOptions_;
            }


        // Evolution

            //- Pre-evolve region
            virtual void preEvolveRegion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

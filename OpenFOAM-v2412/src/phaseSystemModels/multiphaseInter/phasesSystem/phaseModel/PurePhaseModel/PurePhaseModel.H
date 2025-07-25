/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::PurePhaseModel

Description
    Class which represents pure phases, i.e. without any species. Returns an
    empty list of mass fractions.

SourceFiles
    PurePhaseModel.C

\*---------------------------------------------------------------------------*/

#ifndef PurePhaseModel_H
#define PurePhaseModel_H

#include "PtrList.H"
#include "volFields.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiphaseInterSystem;

/*---------------------------------------------------------------------------*\
                           Class PurePhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel, class phaseThermo>
class PurePhaseModel
:
    public BasePhaseModel
{
protected:

    // Protected data

        //- Empty mass fraction field list
        PtrList<volScalarField> Y_;

        //- Empty thermophysical model Ptr
        autoPtr<phaseThermo> thermoPtr_;


public:

    // Constructors

        PurePhaseModel
        (
            const multiphaseInterSystem& fluid,
            const word& phaseName
        );


    //- Destructor
    virtual ~PurePhaseModel() = default;


    // Member Functions

        // Thermo

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& Y();

            //- Solve species fraction equation
            virtual void solveYi
            (
                PtrList<volScalarField::Internal>&,
                PtrList<volScalarField::Internal>&
            );

            //- Access to const thermo
            virtual const phaseThermo& thermo() const;

             //- Access non-const thermo
            virtual phaseThermo& thermo();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PurePhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

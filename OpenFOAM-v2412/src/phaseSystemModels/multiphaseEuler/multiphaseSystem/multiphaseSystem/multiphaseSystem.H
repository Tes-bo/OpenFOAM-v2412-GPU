/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::multiphaseSystem

Description
    Incompressible multi-phase mixture with built in solution for the
    phase fractions with interface compression for interface-capturing.

    Derived from transportModel so that it can be unused in conjunction with
    the incompressible turbulence models.

    Surface tension and contact-angle is handled for the interface
    between each phase-pair.

SourceFiles
    multiphaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_multiphaseEuler_multiphaseSystem_H
#define Foam_multiphaseEuler_multiphaseSystem_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "phaseModel.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dragModel.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class multiphaseSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSystem
:
    public IOdictionary,
    public transportModel
{
public:

    //- Name pair for the interface
    class interfacePair
    :
        public Pair<word>
    {
    public:

        // Ordered hashing (alias)
        using hash = Pair<word>::hasher;

        // Unordered hashing (alias)
        using symmHash = Pair<word>::symmHasher;


        // Constructors

            interfacePair() = default;

            interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
                Pair<word>(alpha1Name, alpha2Name)
            {}

            interfacePair(const phaseModel& alpha1, const phaseModel& alpha2)
            :
                Pair<word>(alpha1.name(), alpha2.name())
            {}


        // Friend Operators

            friend bool operator==
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (0 != Pair<word>::compare(a, b));
            }

            friend bool operator!=
            (
                const interfacePair& a,
                const interfacePair& b
            )
            {
                return (!(a == b));
            }
    };


    typedef
        HashPtrTable
        <
            multiphaseEuler::dragModel, interfacePair, interfacePair::symmHash
        > dragModelTable;

    typedef HashPtrTable<volScalarField, interfacePair, interfacePair::symmHash>
        dragCoeffFields;


private:

    // Private data

        //- Dictionary of phases
        PtrDictionary<phaseModel> phases_;

        const fvMesh& mesh_;
        const surfaceScalarField& phi_;

        volScalarField alphas_;

        typedef HashTable<scalar, interfacePair, interfacePair::symmHash>
            scalarCoeffSymmTable;

        typedef HashTable<scalar, interfacePair, interfacePair::hash>
            scalarCoeffTable;

        scalarCoeffSymmTable sigmas_;
        dimensionSet dimSigma_;

        scalarCoeffSymmTable cAlphas_;

        scalarCoeffTable Cvms_;

        typedef HashTable<dictionary, interfacePair, interfacePair::symmHash>
            interfaceDictTable;

        dragModelTable dragModels_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;


    // Private member functions

        void calcAlphas();

        void solveAlphas();

        tmp<surfaceVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        tmp<surfaceScalarField> nHatf
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;

        void correctContactAngle
        (
            const phaseModel& alpha1,
            const phaseModel& alpha2,
            surfaceVectorField::Boundary& nHatb
        ) const;

        tmp<volScalarField> K
        (
            const phaseModel& alpha1,
            const phaseModel& alpha2
        ) const;


public:

    // Constructors

        //- Construct from components
        multiphaseSystem
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~multiphaseSystem() = default;


    // Member Functions

        //- Return the phases
        const PtrDictionary<phaseModel>& phases() const
        {
            return phases_;
        }

        //- Return the phases
        PtrDictionary<phaseModel>& phases()
        {
            return phases_;
        }

        //- Return the mixture density
        tmp<volScalarField> rho() const;

        //- Return the mixture density for patch
        tmp<scalarField> rho(const label patchi) const;

        //- Return the mixture laminar viscosity
        tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const;

        //- Return the virtual-mass coefficient for the given phase
        tmp<volScalarField> Cvm(const phaseModel& phase) const;

        //- Return the virtual-mass source for the given phase
        tmp<volVectorField> Svm(const phaseModel& phase) const;

        //- Return the table of drag models
        const dragModelTable& dragModels() const
        {
            return dragModels_;
        }

        //- Return the drag coefficients for all of the interfaces
        autoPtr<dragCoeffFields> dragCoeffs() const;

        //- Return the sum of the drag coefficients for the given phase
        tmp<volScalarField> dragCoeff
        (
            const phaseModel& phase,
            const dragCoeffFields& dragCoeffs
        ) const;

        tmp<surfaceScalarField> surfaceTension(const phaseModel& phase) const;

        //- Indicator of the proximity of the interface
        //  Field values are 1 near and 0 away for the interface.
        tmp<volScalarField> nearInterface() const;

        //- Solve for the mixture phase-fractions
        void solve();

        //- Dummy correct
        void correct()
        {}

        //- Read base transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

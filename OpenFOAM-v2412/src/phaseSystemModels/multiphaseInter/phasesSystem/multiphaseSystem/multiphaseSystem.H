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
    Foam::multiphaseSystem

Description
    Class which solves the volume fraction equations for multiple phases

SourceFiles
    multiphaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseSystem_H
#define multiphaseSystem_H

#include "multiphaseInterSystem.H"
#include "UPtrList.H"
#include "phasePairKey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace multiphaseInter
{
/*---------------------------------------------------------------------------*\
                      Class multiphaseSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseSystem
:
    public multiphaseInterSystem
{
protected:

        typedef HashTable<volScalarField::Internal> SuSpTable;

        typedef HashTable<scalar, phasePairKey, phasePairKey::hash>
            scalarTable;

        typedef HashTable<surfaceScalarField> compressionFluxTable;


    // Protected data

        //- Unallocated phase list
        UPtrList<phaseModel> phases_;

        //- Table for compression factors between phases
        scalarTable cAlphas_;

        //- Maximum volumen rate change
        dimensionedScalar ddtAlphaMax_;

        //- Compression fluxed for phases
        compressionFluxTable limitedPhiAlphas_;

        //- Su phase source terms
        SuSpTable Su_;

        //- Sp phase source terms
        SuSpTable Sp_;


    // Protected members

        //- Calculate Sp and Su
        void calculateSuSp();

        //- Solve alphas
        void solveAlphas();


public:

    //- Runtime type information
    TypeName("multiphaseSystem");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            multiphaseSystem,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from fvMesh
        multiphaseSystem(const fvMesh&);


    //- Destructor
    virtual ~multiphaseSystem() = default;


    // Selectors

        static autoPtr<multiphaseSystem> New(const fvMesh& mesh);


     // Member Functions

        //- Solve for the phase fractions
        virtual void solve();


        // Access

            //- Return phases
            const UPtrList<phaseModel>& phases() const;

            //- Return phases
            UPtrList<phaseModel>& phases();

            //- Constant access phase model i
            const phaseModel& phase(const label i) const;

            //- Access phase model i
            phaseModel& phase(const label i);

            //- Access to ddtAlphaMax
            dimensionedScalar ddtAlphaMax() const;

            //- Maximum diffusion number
            scalar maxDiffNo() const;

            //- Access to compression fluxes for phaes
            const compressionFluxTable& limitedPhiAlphas() const;

            //- Access Su
            SuSpTable& Su();

            //- Access Sp
            SuSpTable& Sp();

            //- Read thermophysical properties dictionary
            virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace multiphaseInter

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

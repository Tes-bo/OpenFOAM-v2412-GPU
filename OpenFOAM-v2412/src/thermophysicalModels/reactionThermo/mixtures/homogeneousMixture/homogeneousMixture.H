/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020      OpenCFD Ltd.
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
    Foam::homogeneousMixture

Group
    grpReactionThermophysicalMixtures

Description
    The homogeneous mixture contains species ("b").

SourceFiles
    homogeneousMixture.C

\*---------------------------------------------------------------------------*/

#ifndef homogeneousMixture_H
#define homogeneousMixture_H

#include "basicCombustionMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class homogeneousMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class homogeneousMixture
:
    public basicCombustionMixture
{
    // Private Data

        ThermoType reactants_;
        ThermoType products_;

        mutable ThermoType mixture_;

        //- Regress variable
        volScalarField& b_;

        //- No copy construct
        homogeneousMixture(const homogeneousMixture<ThermoType>&) = delete;


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        homogeneousMixture
        (
            const dictionary& thermoDict,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~homogeneousMixture() = default;


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "homogeneousMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& mixture(const scalar) const;

        const ThermoType& cellMixture(const label celli) const
        {
            return mixture(b_[celli]);
        }

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const
        {
            return mixture(b_[celli]);
        }

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const
        {
            return mixture(b_.boundaryField()[patchi][facei]);
        }


        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const
        {
            return mixture(b_.boundaryField()[patchi][facei]);
        }


        const ThermoType& cellReactants(const label) const
        {
            return reactants_;
        }

        const ThermoType& patchFaceReactants(const label, const label) const
        {
            return reactants_;
        }

        const ThermoType& cellProducts(const label) const
        {
            return products_;
        }

        const ThermoType& patchFaceProducts(const label, const label) const
        {
            return products_;
        }

        //- Read dictionary
        void read(const dictionary&);

        //- Return thermo based on index
        const ThermoType& getLocalThermo(const label speciei) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "homogeneousMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

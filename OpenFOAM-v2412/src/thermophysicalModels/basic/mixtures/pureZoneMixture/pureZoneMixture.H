/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021,2022 OpenCFD Ltd.
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
    Foam::pureZoneMixture

Description
    Version of pureMixture that uses different mixtures for different
    cellZones.
    Every cellZone has to provide an entry for its mixture. A reserved
    entry 'none' is for all unzoned cells.

    Example of the zone based mixture specification:
    \verbatim
    mixture
    {
        solid1
        {
            specie
            {
                molWeight   50;
            }

            transport
            {
                kappa   80;
            }

            thermodynamics
            {
                Hf      0;
                Cp      450;
            }

            equationOfState
            {
                rho     8000;
            }
        }
        solid2
        {
            //- Start off from 'solid1' properties
            ${solid1}

            //- Selectively overwrite properties
            transport
            {
                kappa   8;
            }
        }
    }
    \endverbatim


SourceFiles
    pureZoneMixture.C

\*---------------------------------------------------------------------------*/

#ifndef pureZoneMixture_H
#define pureZoneMixture_H

#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class pureZoneMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class pureZoneMixture
:
    public basicMixture
{
    // Private data

        const fvMesh& mesh_;

        //- Inverse zone info
        labelList zoneID_;

        //- Species data
        PtrList<ThermoType> speciesData_;

        //- Temporary storage for the cell/face mixture thermo data
        mutable ThermoType mixture_;


    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const ThermoType& constructSpeciesData(const dictionary& thermoDict);

        //- No copy construct
        pureZoneMixture(const pureZoneMixture<ThermoType>&) = delete;


public:

    //- The type of thermodynamics this mixture is instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct from dictionary, mesh and phase name
        pureZoneMixture
        (
            const dictionary& thermoDict,
            const fvMesh& mesh,
            const word& phaseName
        );


    //- Destructor
    virtual ~pureZoneMixture() = default;


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "pureZoneMixture<" + ThermoType::typeName() + '>';
        }

        const ThermoType& cellMixture(const label celli) const;

        const ThermoType& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        const ThermoType& cellVolMixture
        (
            const scalar p,
            const scalar T,
            const label celli
        ) const;

        const ThermoType& patchFaceVolMixture
        (
            const scalar p,
            const scalar T,
            const label patchi,
            const label facei
        ) const;

        //- Read dictionary
        void read(const dictionary&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "pureZoneMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::waveModels::cnoidal

Description
    Cnoidal wave model

\*---------------------------------------------------------------------------*/

#ifndef waveModels_cnoidal_H
#define waveModels_cnoidal_H

#include "regularWaveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                           Class cnoidal Declaration
\*---------------------------------------------------------------------------*/

class cnoidal
:
    public regularWaveModel
{
private:

    // Private Member Functions

        //- Initialise
        void initialise
        (
            const scalar H,
            const scalar d,
            const scalar T,
            scalar& mOut,
            scalar& LOut
        ) const;

        //- Wave height
        scalar eta
        (
            const scalar H,
            const scalar Kx,
            const scalar x,
            const scalar Ky,
            const scalar y,
            const scalar omega,
            const scalar t,
            const scalar phase
        ) const;

        scalar eta1D
        (
            const scalar H,
            const scalar m,
            const scalar t,
            const scalar T
        ) const;

        scalar etaMeanSq
        (
            const scalar H,
            const scalar m,
            const scalar T
        ) const;

        vector dEtaDx
        (
            const scalar H,
            const scalar m,
            const scalar uCnoidal,
            const scalar L,
            const scalar K,
            const scalar E
        ) const;

        //- Wave velocity
        vector Uf
        (
            const scalar H,
            const scalar h,
            const scalar Kx,
            const scalar x,
            const scalar Ky,
            const scalar y,
            const scalar omega,
            const scalar t,
            const scalar phase,
            const scalar z
        ) const;


protected:

    // Protected data

        //- `m' coefficient
        scalar m_;


    // Protected Member Functions


        //- Set the water level
        virtual void setLevel
        (
            const scalar t,
            const scalar tCoeff,
            scalarField& level
        ) const;

        //- Calculate the wave model velocity
        virtual void setVelocity
        (
            const scalar t,
            const scalar tCoeff,
            const scalarField& level
        );


public:

    //- Runtime type information
    TypeName("cnoidal");

    //- Constructor
    cnoidal
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~cnoidal() = default;


    // Public Member Functions

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Info
        virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

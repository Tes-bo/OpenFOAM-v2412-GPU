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
    Foam::waveModels::StokesV

Description
    Stokes V wave model

    Reference
    \verbatim
        Skjelbreia, L., Hendrickson, J. (1960).
        Fifth order gravity wave theory.
        Proc. Coastal Engineering Conf., pp. 184-196.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef waveModels_StokesV_H
#define waveModels_StokesV_H

#include "StokesIWaveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                          Class StokesV Declaration
\*---------------------------------------------------------------------------*/

class StokesV
:
    public StokesI
{
private:

    // Private Member Functions

        // Model coefficients

            scalar A11(const scalar h, const scalar k) const;

            scalar A13(const scalar h, const scalar k) const;

            scalar A15(const scalar h, const scalar k) const;

            scalar A22(const scalar h, const scalar k) const;

            scalar A24(const scalar h, const scalar k) const;

            scalar A33(const scalar h, const scalar k) const;

            scalar A35(const scalar h, const scalar k) const;

            scalar A44(const scalar h, const scalar k) const;

            scalar A55(const scalar h, const scalar k) const;

            scalar B22(const scalar h, const scalar k) const;

            scalar B24(const scalar h, const scalar k) const;

            scalar B33(const scalar h, const scalar k) const;

            scalar B33k(const scalar h, const scalar k) const;

            scalar B35(const scalar h, const scalar k) const;

            scalar B35k(const scalar h, const scalar k) const;

            scalar B44(const scalar h, const scalar k) const;

            scalar B55(const scalar h, const scalar k) const;

            scalar B55k(const scalar h, const scalar k) const;

            scalar C1(const scalar h, const scalar k) const;

            scalar C1k(const scalar h, const scalar k) const;

            scalar C2(const scalar h, const scalar k) const;

            scalar C2k(const scalar h, const scalar k) const;

            scalar C3(const scalar h, const scalar k) const;

            scalar C4(const scalar h, const scalar k) const;


        //- Model initialisation
        void initialise
        (
            const scalar H,
            const scalar d,
            const scalar T,
            scalar& kOut,
            scalar& LambdaOut,
            scalar& f1Out,
            scalar& f2Out
        ) const;

        //- Wave height
        scalar eta
        (
            const scalar h,
            const scalar kx,
            const scalar ky,
            const scalar lambda,
            const scalar T,
            const scalar x,
            const scalar y,
            const scalar t,
            const scalar phase
        ) const;

        //- Wave velocity
        vector Uf
        (
            const scalar d,
            const scalar kx,
            const scalar ky,
            const scalar lambda,
            const scalar T,
            const scalar x,
            const scalar y,
            const scalar t,
            const scalar phase,
            const scalar z
        ) const;


protected:

    // Protected Data

        //-
        scalar lambda_;


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
    TypeName("StokesV");

    //- Constructor
    StokesV
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~StokesV() = default;


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

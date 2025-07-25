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
    Foam::waveModels::irregularMultiDirectional

Description
    Multi-directional irregular wave model

\*---------------------------------------------------------------------------*/

#ifndef waveModels_irregularMultiDirectional_H
#define waveModels_irregularMultiDirectional_H

#include "irregularWaveModel.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                  Class irregularMultiDirectional Declaration
\*---------------------------------------------------------------------------*/

class irregularMultiDirectional
:
    public irregularWaveModel
{
private:

    // Proteced Data

        //- Wave periods for irregularMultiDirectional case [seconds]
        List<scalarList> irregWavePeriods_;

        //- Wave heights for irregularMultiDirectional case [metres]
        List<scalarList> irregWaveHeights_;

        //- Wave lengths for irregularMultiDirectional case [metres]
        List<scalarList> irregWaveLengths_;

        //- Wave phases for irregularMultiDirectional case [radians]
        List<scalarList> irregWavePhases_;

        //- Direction of propagation for irregularMultiDirectional case
        //- [degrees] from X axis
        List<scalarList> irregWaveDirs_;


    // Private Member Functions

        //- First Order Wave height (same as StokesI)
        virtual scalar eta
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


protected:

    // Protected Member Functions

        //- Return the wavelength
        virtual scalar waveLength
        (
            const scalar h,
            const scalar T
        ) const;

        //- Wave velocity
        virtual vector Uf
        (
            const scalar d,
            const scalar x,
            const scalar y,
            const scalar t,
            const scalar z
        ) const;

        //-
        virtual vector uMultiDirec
        (
            const scalar irregH,
            const scalar irregWaveOmega,
            const scalar phaseTot,
            const scalar irregWaveKs,
            const scalar z,
            const scalar h,
            const scalar irregDir
        ) const;

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
    TypeName("irregularMultiDirectional");

    //- Constructor
    irregularMultiDirectional
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~irregularMultiDirectional() = default;


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

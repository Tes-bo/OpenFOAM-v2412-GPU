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
    Foam::waveModels::streamFunction

Description
    streamFunction wave model

\*---------------------------------------------------------------------------*/

#ifndef waveModels_streamFunction_H
#define waveModels_streamFunction_H

#include "regularWaveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                       Class streamFunction Declaration
\*---------------------------------------------------------------------------*/

class streamFunction
:
    public regularWaveModel
{
private:

    // Private Member Functions

        //- Wave height
        virtual scalar eta
        (
            const scalar h,
            const scalar kx,
            const scalar ky,
            const scalar T,
            const scalar x,
            const scalar y,
            const scalar omega,
            const scalar t,
            const scalar phase
        ) const;

        //- Wave velocity
        virtual vector Uf
        (
            const scalar d,
            const scalar kx,
            const scalar ky,
            const scalar T,
            const scalar x,
            const scalar y,
            const scalar omega,
            const scalar t,
            const scalar phase,
            const scalar z
        ) const;


protected:

    // Protected data

        //- Mean fluid speed in frame of reference (stream function)
        scalar uMean_;

        //- Stream Function Bj coefficients
        scalarList Bjs_;

        //- Stream Function Ej coefficients
        scalarList Ejs_;


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
    TypeName("streamFunction");

    //- Constructor
    streamFunction
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~streamFunction() = default;


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

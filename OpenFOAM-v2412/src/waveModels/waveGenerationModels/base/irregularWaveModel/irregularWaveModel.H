/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::waveModels::irregularWaveModel

Description

\*---------------------------------------------------------------------------*/

#ifndef waveModels_irregularWaveModel_H
#define waveModels_irregularWaveModel_H

#include "waveGenerationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                     Class irregularWaveModel Declaration
\*---------------------------------------------------------------------------*/

class irregularWaveModel
:
    public waveGenerationModel
{
        //- No copy construct
        irregularWaveModel(const irregularWaveModel&) = delete;

        //- No copy assignment
        void operator=(const irregularWaveModel&) = delete;


protected:

    // Protected data

        //- Ramp time
        scalar rampTime_;

    // Protected Member Functions

        //- Return the time scaling coefficient
        virtual scalar timeCoeff(const scalar t) const;


public:

    //- Runtime type information
    TypeName("irregularWaveModel");


    //- Constructor
    irregularWaveModel
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~irregularWaveModel() = default;


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

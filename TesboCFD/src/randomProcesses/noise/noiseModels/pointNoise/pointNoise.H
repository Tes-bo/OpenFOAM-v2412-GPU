/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
    Foam::noiseModels::pointNoise

Description
    Perform noise analysis on point-based pressure data.

    Input data is read from a dictionary, e.g.

    \verbatim
    // Pressure reference
    pRef            0;

    // Number of samples in sampling window, default = 2^16 (=65536)
    N               4096;

    // Lower frequency bounds
    fl              25;

    // Upper frequency bounds
    fu              25;

    // Start time
    startTime       0;

    windowModel     <modelType>
    <modelType>Coeffs
    {
        ...
    }

    // Pressure data supplied in CSV file format
    file            "pressureData";
    //files           ("pressureData1" "pressureData2");
    nHeaderLine     1;
    refColumn       0;
    componentColumns (1);
    separator       " ";
    mergeSeparators yes;

    \endverbatim

SourceFiles
    pointNoise.C

SeeAlso
    noiseModel.H
    windowModel.H

\*---------------------------------------------------------------------------*/

#ifndef noiseModels_pointNoise_H
#define noiseModels_pointNoise_H

#include "noiseModel.H"
#include "CSV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace noiseModels
{

/*---------------------------------------------------------------------------*\
                         Class pointNoise Declaration
\*---------------------------------------------------------------------------*/

class pointNoise
:
    public noiseModel
{

protected:

    // Protected data

        //- Input file names - optional
        List<fileName> inputFileNames_;


    // Protected Member Functions

        void filterTimeData
        (
            const scalarField& t0,
            const scalarField& p0,
            scalarField& t,
            scalarField& p
        ) const;

        //- Process the CSV data
        void processData
        (
            const label dataseti,
            const Function1Types::CSV<scalar>& data
        );


public:

    //- Runtime type information
    TypeName("pointNoise");

    //- Constructor
    pointNoise
    (
        const dictionary& dict,
        const objectRegistry& obr,
        const word& name = typeName,
        const bool readFields = true
    );

    //- Destructor
    virtual ~pointNoise() = default;


    // Public Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace noiseModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

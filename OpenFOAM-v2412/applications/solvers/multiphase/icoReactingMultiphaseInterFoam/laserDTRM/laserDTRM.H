/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::radiation::laserDTRM

Group
    grpRadiationModels

Description
    Discrete Tray Radiation Method for collimated radiation flux.
    At the moment the particles are injected on a 2D disk in polar
    coordinates providing starting and finish locations.

    The ray is considered planar (no area divergence is considered). This
    model was developed to represent a collimated laser beam.

    The model is developed to work in a VOF framework, representing different
    phases.

    Reflection models are possible to set up between phases.

    Different options are available in order to specify the power distribution
    at the origin of the laser.

SourceFiles
    laserDTRM.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModellaserDTRM_H
#define radiationModellaserDTRM_H

#include "radiationModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Cloud.H"
#include "DTRMParticle.H"
#include "reflectionModel.H"
#include "Function1.H"
#include "interpolation2DTable.H"
#include "labelField.H"
#include "phasePairKey.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class laserDTRM Declaration
\*---------------------------------------------------------------------------*/

class laserDTRM
:
    public radiationModel
{
public:

    //  Public enumeration

        //- Modes of power distribution
        enum powerDistributionMode
        {
            pdGaussian,
            pdManual,
            pdUniform,
            pdGaussianPeak
        };


private:


    // Private types

        typedef
            HashTable<dictionary, phasePairKey, phasePairKey::hash> dictTable;


        typedef
            HashTable
            <
                autoPtr<reflectionModel>,
                phasePairKey,
                phasePairKey::hash
            > reflectionModelTable;


    // Private data

        static const Enum<powerDistributionMode> powerDistNames_;

        //- Operating mode for power distribution
        powerDistributionMode mode_;

        //- DTRM particle cloud
        Cloud<DTRMParticle> DTRMCloud_;

        //- Number of particles
        label nParticles_;

        //- Number of theta angles
        label ndTheta_;

        //- Number of radials
        label ndr_;

        //- Maximum tracking length for particles
        scalar maxTrackLength_;


        // Laser parameters

            //- Focal laser position
            autoPtr<Function1<point>> focalLaserPosition_;

            //- Direction
            autoPtr<Function1<vector>> laserDirection_;

            //- Focal lase radius
            scalar focalLaserRadius_;

            //- Quality beam laser (tan of the angle between the main direction
            //  (laserDirection) and the spreading rays (0: flat, 1:fully)
            scalar qualityBeamLaser_;


            // Gaussian power distribution entries

                //- Sigma deviation
                scalar sigma_;

                //- Peak flux for GaussianPeak mode
                scalar I0_;

                //- Power
                autoPtr<Function1<scalar>> laserPower_;


            // Manual power distribution entries

                //- Manual power distribution table. (theta, radius)
                autoPtr<interpolation2DTable<scalar>> powerDistribution_;


            // Reflection sub-model

                //- Reflection model
                reflectionModelTable reflections_;

                //- Reflection switch
                bool reflectionSwitch_;

                //- Alpha value at which reflection is set
                scalar alphaCut_;


        // Fields


            //- Absorption coefficient
            volScalarField a_;

            //- Emission coefficient
            volScalarField e_;

            //- Emission contribution
            volScalarField E_;

            //- Volumetric heat source [W/m3]
            volScalarField Q_;


    // Private Member Functions

        //- Initialise
        void initialise();

        //- Initialise reflection model
        void initialiseReflection();

        //- Calculate Intensity of the laser at p(t, theta) [W/m2]
        scalar calculateIp(scalar r, scalar theta);

        // Normal surface on the interface
        tmp<volVectorField> nHatfv
        (
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ) const;


        //- No copy construct
        laserDTRM(const laserDTRM&) = delete;

        //- No copy assignment
        void operator=(const laserDTRM&) = delete;


public:

    //- Runtime type information
    TypeName("laserDTRM");


    // Constructors

        //- Construct from components
        laserDTRM(const volScalarField& T);

        //- Construct from components
        laserDTRM(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~laserDTRM() = default;


    // Member functions

    // Edit

        //- Solve radiation equation(s)
        void calculate();

        //- Read radiation properties dictionary
        bool read();

        //- Number of bands for this radiation model
        virtual label nBands() const;


    // Access

        //- Source term component (for power of T^4)
        virtual tmp<volScalarField> Rp() const;

        //- Source term component (constant)
        virtual tmp<DimensionedField<scalar, volMesh>> Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

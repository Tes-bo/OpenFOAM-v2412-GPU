/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::fv::rotorDiskSource

Group
    grpFvOptionsSources

Description
    Applies cell-based momentum sources on velocity (i.e. \c U)
    within a specified cylindrical region to approximate the mean
    effects of rotor forces.

    Sources applied to:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity    [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    rotorDiskSource1
    {
        // Mandatory entries (unmodifiable)
        type            rotorDiskSource;

        // Mandatory entries (runtime modifiable)
        fields          (\<field1\> \<field2\> ... \<fieldN\>);
        rhoRef          1.0;
        rpm             5.5;
        nBlades         3;
        inletFlowType   local;
        geometryMode    auto;
        tipEffect       0.96; // [0, 1]
        refDirection    (-1 0 0 );
        flapCoeffs
        {
            beta0       0.0;
            beta1c      0.0;
            beta2s      0.0;
        }

        trimModel       fixed;

        blade
        {
            // see bladeModel.H for documentation
        }

        profiles
        {
            profile1
            {
                type    lookup; // lookup || series
                ...
                // see lookupProfile.H or seriesProfile.H for documentation
            }
            profile2
            {
                ...
            }
        }

        // Conditional mandatory entries (runtime modifiable)

            // when inletFlowType=fixed
            inletVelocity    (1 0 0);

            // when inletFlowType=surfaceNormal
            inletNormalVelocity  5.0;

            // when geometryMode=auto
            pointAbove       (1 0 0);

            // when geometryMode=specified
            origin           (0 0 0);
            axis             (1 0 0);

        // Mandatory/Optional (inherited) entries
        ...
    }

    where the entries mean:
    \table
      Property   | Description                          | Type  | Reqd | Dflt
      type       | Type name: rotorDiskSource           | word  | yes  | -
      fields     | Names of operand fields              | wordList | yes | -
      rhoRef     | Reference density for incompressible case | scalar | yes | -
      rpm        | Rotational speed [rad/s]             | scalar | yes | -
      nBlades    | Number of rotor blades               | label | yes  | -
      tipEffect  | Ratio of blade radius beyond which lift=0 | scalar | yes | -
      refDirection | Reference direction used as reference <!--
                 --> for psi angle                      | vector | yes | -
      beta0      | Coning angle [deg]                   | scalar | yes | -
      beta1c  | Lateral flapping coeff (cos coeff) [deg]      | scalar | yes | -
      beta2s  | Longitudinal flapping coeff (sin coeff) [deg] | scalar | yes | -
      inletFlowType | Inlet flow type specification     | word  | yes  | -
      inletVelocity | Inlet velocity [m/s]              | vector | cndtnl | -
      inletNormalVelocity | Inlet normal velocity magnitude <!--
                                                    --> | scalar | cndtnl | -
      geometryMode  | Geometry mode specification       | word  | yes  | -
      pointAbove  | Correct the axis direction using a point <!--
                  --> above the rotor                    | vector | cndtnl | -
      origin | Origin of the specified coordinate system | vector | cndtnl | -
      axis   | Axis of the specified coordinate system   | vector | cndtnl | -
      trimModel  | Trim model specification             | word  | yes  | -
      blade      | Blade properties                | dictionary | yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink
     - \link trimModel.H \endlink
     - \link bladeModel.H \endlink
     - \link profileModel.H \endlink

    Options for the \c inletFlowType entry:
    \verbatim
      fixed    | Use a specified velocity
      local    | Use local flow conditions
      surfaceNormal | Use a specified normal velocity (positive towards rotor)
    \endverbatim

    Options for the \c geometryMode entry:
    \verbatim
      auto      | Determine rotor coordinate system from selected cells
      specified | Use a specified coordinate system
    \endverbatim

    Options for the \c trimModel entry:
    \verbatim
      fixed       | Fixed blade/angle characteristics
      targetForce | Target thrust/torque by using time-variant blade characs.
    \endverbatim

Note
  - The entry \c rpm is positive anti-clockwise
    when looking along -ve lift direction.

See also
  - Foam::trimModel
  - Foam::bladeModel
  - Foam::lookupProfile
  - Foam::seriesProfile

SourceFiles
    rotorDiskSource.C
    rotorDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef rotorDiskSource_H
#define rotorDiskSource_H

#include "cellSetOption.H"
#include "cylindricalCS.H"
#include "Enum.H"
#include "bladeModel.H"
#include "profileModelList.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class trimModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class rotorDiskSource Declaration
\*---------------------------------------------------------------------------*/

class rotorDiskSource
:
    public fv::cellSetOption
{
public:

    // Public Enumerations

        //- Options for the geometry type specification
        enum geometryModeType
        {
            gmAuto,
            gmSpecified
        };

        //- Names for geometryModeType
        static const Enum<geometryModeType> geometryModeTypeNames_;

        //- Options for the inlet flow type specification
        enum inletFlowType
        {
            ifFixed,
            ifSurfaceNormal,
            ifLocal
        };

        //- Names for inletFlowType
        static const Enum<inletFlowType> inletFlowTypeNames_;


protected:

    // Helper structures to encapsulate flap and trim data
    // Note: all input in degrees (converted to radians internally)

        struct flapData
        {
            scalar beta0;   // coning angle
            scalar beta1c;  // lateral flapping coeff (cos coeff)
            scalar beta2s;  // longitudinal flapping coeff (sin coeff)
        };


    // Protected Data

        //- Reference density for incompressible case
        scalar rhoRef_;

        //- Rotational speed [rad/s]
        //  Positive anti-clockwise when looking along -ve lift direction
        scalar omega_;

        //- Number of blades
        label nBlades_;

        //- Inlet flow type
        inletFlowType inletFlow_;

        //- Inlet velocity for specified inflow
        vector inletVelocity_;

        //- Tip effect [0-1]
        //  Ratio of blade radius beyond which lift=0
        scalar tipEffect_;

        //- Blade flap coefficients [rad/s]
        flapData flap_;

        //- Cell centre positions in local rotor frame
        //  (Cylindrical r-theta-z)
        List<point> x_;

        //- Rotation tensor for flap angle
        List<tensor> Rcone_;

        //- Area [m2]
        List<scalar> area_;

        //- Rotor local cylindrical coordinate system (r-theta-z)
        coordSystem::cylindrical coordSys_;

        //- Cached rotation tensors for cylindrical coordinates
        autoPtr<tensorField> Rcyl_;

        //- Maximum radius
        scalar rMax_;

        //- Trim model
        autoPtr<trimModel> trim_;

        //- Blade data
        bladeModel blade_;

        //- Profile data
        profileModelList profiles_;


    // Protected Member Functions

        //- Check data
        void checkData();

        //- Set the face areas per cell, and optionally correct the rotor axis
        void setFaceArea(vector& axis, const bool correct);

        //- Create the coordinate system
        void createCoordinateSystem();

        //- Construct geometry
        void constructGeometry();

        //- Return the inlet flow field
        tmp<vectorField> inflowVelocity(const volVectorField& U) const;

        //- Helper function to write rotor values
        template<class Type>
        void writeField
        (
            const word& name,
            const List<Type>& values,
            const bool writeNow = false
        ) const;


public:

    //- Runtime type information
    TypeName("rotorDisk");


    // Constructors

        //- Construct from components
        rotorDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        rotorDiskSource(const rotorDiskSource&) = delete;

        //- No copy assignment
        void operator=(const rotorDiskSource&) = delete;


    //- Destructor
    virtual ~rotorDiskSource() = default;


    // Member Functions

        // Access

            //- Return the reference density for incompressible case
            inline scalar rhoRef() const;

            //- Return the rotational speed [rad/s]
            //  Positive anti-clockwise when looking along -ve lift direction
            inline scalar omega() const;

            //- Return the cell centre positions in local rotor frame
            //  (Cylindrical r-theta-z)
            inline const List<point>& x() const;

            //- Return the rotor coordinate system (r-theta-z)
            inline const coordSystem::cylindrical& coordSys() const;


        // Evaluation

            //- Calculate forces
            template<class RhoFieldType>
            void calculate
            (
                const RhoFieldType& rho,
                const vectorField& U,
                const scalarField& thetag,
                vectorField& force,
                const bool divideVolume = true,
                const bool output = true
            ) const;


        // Source term addition

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rotorDiskSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rotorDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

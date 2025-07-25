/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::functionObjects::setFlow

Group
    grpFieldFunctionObjects

Description
    Provides options to set the velocity and flux fields as a function of time.

    Useful for testing advection behaviour of numerical schemes by e.g.
    imposing solid body rotation, vortex flows.  All types include a scaling
    \c Function1 type enabling the strength of the transformation to vary
    as a function of time.

    Operands:
    \table
      Operand       | Type                       | Location
      input         | {vol,surface}\<Type\>Field <!--
                --> | $FOAM_CASE/\<time\>/\<inpField\>
      output file   | -                          | -
      output field  | {vol,surface}\<Type\>Field <!--
                --> | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

    where \c \<Type\>=Scalar/Vector/SphericalTensor/SymmTensor/Tensor.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    setFlow1
    {
        // Mandatory entries (unmodifiable)
        type        setFlow;
        libs        (fieldFunctionObjects);

        // Mandatory entries (runtime modifiable)
        mode        rotation;
        scale       1;

        // Optional entries (runtime modifiable)
        U           U;
        rho         none;
        phi         phi;
        reverseTime 1;

        // When mode=function
        velocity    (1 0 0);

        // When mode=rotation
        omega       6.28318530718;
        origin      (0.5 0 0.5);
        refDir      (1 0 0);
        axis        (0 1 0);

        // When mode=(vortex2D|vortex3D)
        origin      (0.5 0 0.5);
        refDir      (1 0 0);
        axis        (0 1 0);

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                         | Type | Req'd | Dflt
      type         | Type name: readFields               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects  | word |  yes  | -
      mode         | Operating mode - see below          | word |  yes  | -
      scale        | Scaling function        | Function1<scalar> | yes  | -
      U            | Name of velocity field              | word |  no   | U
      rho          | Name of density field               | word |  no   | none
      phi          | Name of flux field                  | word |  no   | phi
      reverseTime  <!--
      --> | Specified time to reverse flow direction     | scalar  | no | VGREAT
      velocity     <!--
      --> | Velocity function           | Function1<vector> | conditional | -
      omega        <!--
      --> | Rotational speed function   | Function1<scalar> | conditional | -
      origin       <!--
      --> | Rotation vector origin      | vector            | conditional | -
      refDir       <!--
      --> | Rotation vector reference direction   | vector  | conditional | -
      axis         <!--
      --> | Rotation vector axis                  | vector  | conditional | -
    \endtable

    Options for the \c mode entry:
    \verbatim
      function
      rotation
      vortex2D
      vortex3D
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - ExtendedCodeGuide::functionObjects::field::setFlow

SourceFiles
    setFlow.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_setFlow_H
#define functionObjects_setFlow_H

#include "fvMeshFunctionObject.H"
#include "Function1.H"
#include "Enum.H"
#include "point.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class setFlow Declaration
\*---------------------------------------------------------------------------*/

class setFlow
:
    public fvMeshFunctionObject
{
    // Private Enumerations

        //- Options for the operating mode
        enum class modeType
        {
            FUNCTION,
            ROTATION,
            VORTEX2D,
            VORTEX3D
        };

        //- Names for modeType
        static const Enum<modeType> modeTypeNames;


   // Private Data

        //- Operating mode
        modeType mode_;

        //- Name of velocity field
        word UName_;

        //- Name of density field
        word rhoName_;

        //- Name of flux field
        word phiName_;

        //- Reverse time
        scalar reverseTime_;

        //- Scaling function
        autoPtr<Function1<scalar>> scalePtr_;


        // Rotation

            //- Origin
            point origin_;

            //- Rotation tensor for rotational mode
            tensor R_;

            //- Rotational speed function
            autoPtr<Function1<scalar>> omegaPtr_;


        // Function

            //- Velocity function
            autoPtr<Function1<vector>> velocityPtr_;


    // Private Member Functions

        //- Set the flux field
        void setPhi(const volVectorField& U);


public:

    //- Runtime type information
    TypeName("setFlow");


    // Constructors

        //- Construct from Time and dictionary
        setFlow
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        setFlow(const setFlow&) = delete;

        //- No copy assignment
        void operator=(const setFlow&) = delete;


    //- Destructor
    virtual ~setFlow() = default;


    // Member Functions

        //- Read the setFlow data
        virtual bool read(const dictionary& dict);

        //- Do nothing
        virtual bool execute();

        //- Calculate the setFlow and write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

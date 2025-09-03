/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::functionObjects::stabilityBlendingFactor

Group
    grpFieldFunctionObjects

Description
    Computes the \c stabilityBlendingFactor to be used by the
    local blended convection scheme. The output is a surface field weight
    between 0-1.

    The weight of a blended scheme, i.e. \c w, is given by a function of
    the blending factor, \c f:

    \f[
        w = f_{scheme_1} + (1 - f_{scheme_2})
    \f]

    The factor is calculated based on six criteria:
    \verbatim
      1. mesh non-orthogonality field
      2. magnitude of cell centres gradient
      3. convergence rate of residuals
      4. faceWeight
      5. skewness
      6. Courant number
    \endverbatim

    The user can enable them individually.

    For option 1, the following relation is used, where \f$\phi_1\f$ is
    the non-orthogonality:
    \f[
        fNon =
            min
            (
                max
                (
                    0.0,
                    (\phi_1 - max(\phi_1))
                    /(min(\phi_1) - max(\phi_1))
                ),
                1.0
            )
    \f]

    For option 2, the following relation is used, where \f$\phi_2\f$ is
    the magnitude of cell centres gradient (Note that \f$\phi_2 = 3\f$
    for orthogonal meshes):

    \f[
        fMagGradCc =
            min
            (
                max
                (
                    0.0,
                    (\phi_2 - max(\phi_2))
                    / (min(\phi_2) - max(\phi_2))
                ),
                1.0
            )
    \f]

    For option 3, a PID control is used in order to control residual
    unbounded fluctuations for individual cells.

    \f[
        factor =
            P*residual
            + I*residualIntegral
            + D*residualDifferential
    \f]

    where \c P, \c I and \c D are user inputs.

    The following relation is used:
    \f[
        fRes = (factor - meanRes)/(maxRes*meanRes);
    \f]

    where
    \vartable
        meanRes | Average(residual)
        maxRes  | User input
    \endvartable

    Note that \f$f_{Res}\f$ will blend more towards one as
    the cell residual is larger then the domain mean residuals.


    For option 4, the following relation is used, where \f$\phi_4\f$ is
    the face weight (Note that \f$\phi_4 = 0.5\f$ for orthogonal meshes):

    \f[
        ffaceWeight = min
        (
            max
            (
                0.0,
                (min(\phi_4) - \phi_4)
                / (min(\phi_4) - max(\phi_4))
            ),
            1.0
        )
    \f]


    For option 5, the following relation is used, where \f$\phi_5\f$ is
    the cell skewness:

    \f[
        fskewness =
        min
        (
            max
            (
                0.0,
                (\phi_5    - max(\phi_5))
                / (min(\phi_5) - max(\phi_5))
            ),
            1.0
        )
    \f]


    For option 6, the following relation is used:

    \f[
        fCoWeight = clamp((Co - Co1)/(Co2 - Co1), zero_one{});
    \f]

    where
    \vartable
        Co1 | Courant number below which scheme2 is used
        Co2 | Courant number above which scheme1 is used
    \endvartable

    The final factor is determined by:

    \f[
        f = max(fNon, fMagGradCc, fRes, ffaceWeight, fskewness, fCoWeight)
    \f]

    An indicator (volume) field, named \c blendedIndicator
    is generated if the log flag is on:
    - 1 represent scheme1 as active,
    - 0 represent scheme2 as active.

    Additional reporting is written to the standard output, providing
    statistics as to the number of cells used by each scheme.

    Operands:
    \table
      Operand        | Type           | Location
      input          | -              | -
      output file    | dat | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field   | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    stabilityBlendingFactor1
    {
        // Mandatory entries (unmodifiable)
        type                stabilityBlendingFactor;
        libs                (fieldFunctionObjects);

        // Mandatory entries (unmodifiable)
        field               <field>;    // U;
        result              <outField>; // UBlendingFactor;

        // Optional entries (runtime modifiable)
        tolerance           0.001;

        // Any of the options can be chosen in combinations

        // Option-1
        switchNonOrtho      true;
        nonOrthogonality    nonOrthoAngle;
        maxNonOrthogonality 20;
        minNonOrthogonality 60;

        // Option-2
        switchGradCc        true;
        maxGradCc           3;
        minGradCc           4;

        // Option-3
        switchResiduals     true;
        maxResidual         10;
        residual            initialResidual:p;
        P                   1.5;
        I                   0;
        D                   0.5;

        // Option-4
        switchFaceWeight    true;
        maxFaceWeight       0.3;
        minFaceWeight       0.2;

        // Option-5
        switchSkewness      true;
        maxSkewness         2;
        minSkewness         3;

        // Option-6
        switchCo            true;
        U                   U;
        Co1                 1;
        Co2                 10;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    Example of function object specification to calculate the \c residuals used
    by \c stabilityBlendingFactor. The following writes 'initialResidual:p'
    field
    \verbatim
    residuals
    {
        type            residuals;
        libs            (utilityFunctionObjects);
        writeFields     true;
        writeControl    writeTime;
        fields          (p);
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Req'd | Dflt
      type         | Type name: stabilityBlendingFactor    | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      field        | Name of operand field                 | word |  yes  | -
      result  | Name of surface field to be used in the localBlended scheme <!--
          --> | word | yes
      switchNonOrtho | Select non-orthogonal method        | bool | no | false
      nonOrthogonality | Name of the non-orthogonal field  <!--
                   --> | word | no | nonOrthoAngle
      maxNonOrthogonality| Maximum non-orthogonal for scheme2 | scalar | no | 20
      minNonOrthogonality| Minimum non-orthogonal for scheme1 | scalar | no | 60
      switchGradCc | Select cell centre gradient method    | bool | no | false
      maxGradCc| Maximum gradient for scheme2              | scalar | no | 2
      minGradCc| Minimum gradient for scheme1              | scalar | no | 4
      switchResiduals | Select residual evolution method   | bool | no | false
      residual    | Name of the residual field | word | no | initialResidual:p
      maxResidual| Maximum residual-mean ratio for scheme1 | scalar | no | 10
      P       | Proportional factor for PID                | scalar | no | 3
      I       | Integral factor for PID                    | scalar | no | 0
      D       | Differential factor for PID                | scalar | no | 0.25
      switchFaceWeight | Select face weight method         | bool | no | false
      faceWeight | Name of the faceWeight field       | word | no | faceWeight
      maxFaceWeight | Maximum face weight for scheme1      | scalar | no | 0.2
      minFaceWeight | Minimum face weight for scheme2      | scalar | no | 0.3
      switchSkewness   | Select skewness method            | bool | no | false
      skewness | Name of the skewness field           | word | no | skewness
      maxSkewness | Maximum skewness for scheme2           | scalar | no | 2
      minSkewness | Minimum skewness for scheme1           | scalar | no | 3
      switchCo         | Select Co blended method          | bool | no | false
      U   | Name of the flux field for Co blended          | word | no | U
      Co1 | Courant number below which scheme2 is used     | scalar | no | 1
      Co2 | Courant number above which scheme1 is used     | scalar | no | 10
      tolerance    | Tolerance for number of blended cells | scalar | no | 0.001
    \endtable

    The \c result entry is the field which is read by the \c localBlended scheme
    specified in \c fvSchemes. This name is determined by the \c localBlended
    class.

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::writeFile
    - ExtendedCodeGuide::functionObjects::field::stabilityBlendingFactor

SourceFiles
    stabilityBlendingFactor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stabilityBlendingFactor_H
#define functionObjects_stabilityBlendingFactor_H

#include "fieldExpression.H"
#include "writeFile.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class stabilityBlendingFactor Declaration
\*---------------------------------------------------------------------------*/

class stabilityBlendingFactor
:
    public fieldExpression,
    public writeFile
{
    // Private Member Data

        // Switches

            //- Switch for non-orthogonality
            Switch nonOrthogonality_;

            //- Switch for grad of cell centres
            Switch gradCc_;

            //- Switch for residuals
            Switch residuals_;

            //- Switch for face weight
            Switch faceWeight_;

            //- Switch for skewness
            Switch skewness_;

            //- Switch for Co
            Switch Co_;


        // Lower and upper limits

            //- Maximum non-orthogonality for fully scheme 2
            scalar maxNonOrthogonality_;

            //- Minimum non-orthogonality for fully scheme 1
            scalar minNonOrthogonality_;

            //- Maximum  gradcc for fully scheme 2
            scalar maxGradCc_;

            //- Minimum  gradcc for fully scheme 1
            scalar minGradCc_;

            //- Maximum ratio to average residual for scheme 2
            scalar maxResidual_;

            //- Minimum face weight for fully scheme 2
            scalar minFaceWeight_;

            //- Maximum face weight for fully scheme 1
            scalar maxFaceWeight_;

            //- Maximum skewness for fully scheme 2
            scalar maxSkewness_;

            //- Minimum skewness for fully scheme 1
            scalar minSkewness_;

            //- Maximum Co for fully scheme 2
            scalar Co1_;

            //- Minimum Co for fully scheme 1
            scalar Co2_;


        // File names

            //- Name of the non-orthogonality field
            word nonOrthogonalityName_;

            //- Name of the face weight field
            word faceWeightName_;

            //- Name of the skewnes field
            word skewnessName_;

            //- Name of the residual field
            word residualName_;

            //- Name of the U used for Co based blended
            word UName_;


        //- Tolerance used when calculating the number of blended cells
        scalar tolerance_;


        //- Error fields
        scalarField error_;
        scalarField errorIntegral_;
        scalarField oldError_;
        scalarField oldErrorIntegral_;

        //- Proportional gain
        scalar P_;

        //- Integral gain
        scalar I_;

        //- Derivative gain
        scalar D_;


    // Private Member Functions

        //- Init fields
        bool init(bool first);

        //- Return access to the indicator field
        volScalarField& indicator();

        //- Calculate statistics
        void calcStats(label&, label&, label&) const ;

        //- Calculate the blending factor field and return true if successful
        virtual bool calc();


protected:

    // Protected Member Functions

        //- Write the file header
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("stabilityBlendingFactor");


    // Constructors

        //- Construct from Time and dictionary
        stabilityBlendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        stabilityBlendingFactor(const stabilityBlendingFactor&) = delete;

        //- No copy assignment
        void operator=(const stabilityBlendingFactor&) = delete;


    //- Destructor
    virtual ~stabilityBlendingFactor() = default;


    // Member Functions

        //- Read the stabilityBlendingFactor data
        virtual bool read(const dictionary&);

        //- Write the stabilityBlendingFactor
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

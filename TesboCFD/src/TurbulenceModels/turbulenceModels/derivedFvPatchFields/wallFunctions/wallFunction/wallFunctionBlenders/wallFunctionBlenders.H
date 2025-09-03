/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::wallFunctionBlenders

Description
    The class \c wallFunctionBlenders is a base class that
    hosts common entries for various derived wall-function boundary
    conditions to be used in low- and high-Reynolds number applications.

    Reference:
    \verbatim
        Default model coefficients (tag:VM):
            Versteeg, H. K., & Malalasekera, W. (2011).
            An introduction to computational fluid dynamics: The finite
            volume method. Harlow: Pearson Education.
            Subsection "3.5.2 k-epsilon model".

        Binomial blending of the viscous and inertial sublayers (tag:ME):
            Menter, F., & Esch, T. (2001).
            Elements of industrial heat transfer prediction.
            In Proceedings of the 16th Brazilian Congress of Mechanical
            Engineering (COBEM), November 2001. vol. 20, p. 117-127.

        Exponential/Max blending of the viscous and inertial sublayers (tag:PH):
            Popovac, M., & Hanjalić, K. (2007).
            Compound wall treatment for RANS computation of complex
            turbulent flows and heat transfer.
            Flow, turbulence and combustion, 78(2), 177-202.
            DOI:10.1007/s10494-006-9067-x

        Tanh blending of the viscous and inertial sublayers (tag:KAS):
            Knopp, T., Alrutz, T., & Schwamborn, D. (2006).
            A grid and flow adaptive wall-function method for RANS
            turbulence modelling.
            Journal of Computational Physics, 220(1), 19-40.
            DOI:10.1016/j.jcp.2006.05.003
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Top-level entries
        ...

        // Optional entries
        blending        stepwise;
        n               4.0;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type   | Reqd  | Deflt
      blending  | Viscous/inertial sublayer blending | word | no  | stepwise
      n         | Binomial blending exponent     | scalar | no    | 2.0
    \endtable

    Options for the \c blending entry:
    \verbatim
      stepwise    | Stepwise switch (discontinuous)
      max         | Maximum value switch (discontinuous)
      binomial    | Binomial blending (smooth)
      exponential | Exponential blending (smooth)
      tanh        | Tanh blending (smooth)
    \endverbatim

    wherein predictions of a given quantity (e.g. \c nut) for viscous and
    inertial sublayers are blended according to the following expressions:

    - \c stepwise (default):

    \f[
        \phi = \phi_{log} \qquad if \quad y^+ > y^+_{intersection}
    \f]

    \f[
        \phi = \phi_{vis} \qquad if \quad y^+ <= y^+_{intersection}
    \f]

    where
    \vartable
      \phi_{vis} | \f$\phi\f$ prediction in viscous sublayer
      \phi_{log} | \f$\phi\f$ prediction in inertial sublayer
      y^+   | estimated wall-normal height of cell centre in wall units
      y^+_{intersection} | estimated \f$y^+\f$ where sublayers intersect
    \endvartable


    - \c max (PH:Eq. 27):

    \f[
        \nu_t = max({\nu_t}_{vis}, {\nu_t}_{log})
    \f]


    - \c binomial (ME:Eqs. 15-16):

    \f[
        \phi = ((\phi_{vis})^n + (\phi_{log})^n)^{1/n}
    \f]
    where
    \vartable
      n               | Binomial blending exponent
    \endvartable


    - \c exponential (PH:Eq. 32):

    \f[
        \phi = \phi_{vis} \exp[-\Gamma] + \phi_{log} \exp[-1/\Gamma]
    \f]

    where (PH:Eq. 31)
    \vartable
      \Gamma       | Blending expression
      \Gamma       | \f$0.01 (y^+)^4 / (1.0 + 5.0 y^+)\f$
    \endvartable


    - \c tanh (KAS:Eqs. 33-34):

    \f[
        \phi = \eta \phi_{b1} + (1 - \eta)\phi_{b2}
    \f]

    where
    \vartable
        \eta        | \f$tanh((y^+/10)^4)\f$
        \phi_{b1}   | \f$\phi_{vis} + \phi_{log}\f$
        \phi_{b2}   | \f$(\phi_{vis}^{1.2} + \phi_{log}^1.2)^{1/1.2}\f$
    \endvartable

    \c G predictions for viscous and inertial sublayers are blended
    in a stepwise manner, and \c G below \f$y^+_{intersection}\f$
    (i.e. in the viscous sublayer) is presumed to be zero.

    Blending treatments are enabled for the following wall functions:
      - \link epsilonWallFunctionFvPatchScalarField.H \endlink
      - \link omegaWallFunctionFvPatchScalarField.H \endlink
      - \link nutkWallFunctionFvPatchScalarField.H \endlink
      - \link nutUWallFunctionFvPatchScalarField.H \endlink

SourceFiles
    wallFunctionBlenders.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallFunctionBlenders_H
#define Foam_wallFunctionBlenders_H

#include "scalarFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class Ostream;
template<class EnumType> class Enum;

/*---------------------------------------------------------------------------*\
                    Class wallFunctionBlenders Declaration
\*---------------------------------------------------------------------------*/

class wallFunctionBlenders
{
protected:

    // Protected Enumerations

        //- Options for the blending treatment of viscous and inertial sublayers
        enum blenderType
        {
            STEPWISE,       //!< "Stepwise switch (discontinuous)"
            MAX,            //!< "Maximum value switch (discontinuous)"
            BINOMIAL,       //!< "Binomial blending (smooth)"
            EXPONENTIAL,    //!< "Exponential blending (smooth)"
            TANH            //!< "Tanh blending (smooth)"
        };

        //- Names for blenderType
        static const Enum<blenderType> blenderTypeNames;


    // Protected Data

        //- Blending treatment
        enum blenderType blender_;

        //- Binomial blending exponent being used when
        //- blenderType is blenderType::BINOMIAL
        scalar n_;


public:

    // Constructors

        //- Default construct with default coefficients
        wallFunctionBlenders();

        //- Construct from dictionary and default parameters
        wallFunctionBlenders
        (
            const dictionary& dict,
            const blenderType blender,
            const scalar n
        );


    // Member Functions

    // I-O

        //- Write wall-function blending data as dictionary entries
        void writeEntries(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

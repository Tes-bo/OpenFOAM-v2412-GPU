/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::turbulentDigitalFilterInletFvPatchField

Group
    grpInletBoundaryConditions

Description
    Digital-filter based boundary condition for vector- and scalar-based
    quantities (e.g. \c U or \c T) to generate synthetic turbulence-alike
    time-series from input turbulence statistics for LES and DES
    turbulent flow computations.

    References:
    \verbatim
        Digital-filter method-based generator (DFM) (tag:KSJ):
            Klein, M., Sadiki, A., & Janicka, J. (2003).
            A digital filter based generation of inflow data for spatially
            developing direct numerical or large eddy simulations.
            Journal of computational Physics, 186(2), 652-665.
            DOI:10.1016/S0021-9991(03)00090-1

        Forward-stepwise method-based generator (FSM) (tag:XC)
            Xie, Z. T., & Castro, I. P. (2008).
            Efficient generation of inflow conditions for
            large eddy simulation of street-scale flows.
            Flow, turbulence and combustion, 81(3), 449-470.
            DOI:10.1007/s10494-008-9151-5

        Mass-inflow rate correction (tag:KCX):
            Kim, Y., Castro, I. P., & Xie, Z. T. (2013).
            Divergence-free turbulence inflow conditions for
            large-eddy simulations with incompressible flow solvers.
            Computers & Fluids, 84, 56-68.
            DOI:10.1016/j.compfluid.2013.06.001

        Scalar fluctuations (tag:XHW):
            Xie, Z. T., Hayden, P., & Wood, C. R. (2013).
            Large-eddy simulation of approaching-flow stratification
            on dispersion over arrays of buildings.
            Atmospheric Environment, 71, 64-74.
            DOI:10.1016/j.atmosenv.2013.01.054

            Okaze, T., & Mochida, A. (2017).
            Cholesky decomposition–based generation of artificial
            inflow turbulence including scalar fluctuation.
            Computers & Fluids, 159, 23-32.
            DOI:10.1016/j.compfluid.2017.09.005
    \endverbatim

    In \c DFM or \c FSM, a random number set (mostly white noise), and a group
    of target statistics (mostly mean flow, Reynolds-stress tensor profiles and
    integral-scale sets) are merged into a new number set (stochastic
    time-series, yet consisting of the statistics) by a chain of mathematical
    operations whose characteristics are designated by the target statistics,
    so that the realised statistics of the new sets could match the target.

    \verbatim
    Random number sets ---->-|
                             |
                         DFM or FSM ---> New stochastic time-series consisting
                             |           turbulence statistics
    Turbulence statistics ->-|
    \endverbatim

    The main difference between \c DFM and \c FSM is that \c FSM replaces
    the expensive-to-run streamwise convolution summation in \c DFM by a simpler
    and an almost-equivalent-in-effect numerical procedure in order to reduce
    computational costs. Accordingly, \c FSM potentially brings computational
    resource advantages for computations involving relatively large streamwise
    integral-scale sets and small time-steps.

    Synthetic turbulence is generated on a virtual rectangular structured-mesh
    plane and is mapped onto this patch by the selected mapping method.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type        turbulentDigitalFilterInlet;
        n           (<label> <label>);
        mean        <PatchFunction1<vector>> or <PatchFunction1<scalar>>;
        R           <PatchFunction1<symmTensor>> or <PatchFunction1<scalar>>;
        L           <tensor> or <vector>;

        // Optional entries
        kernelType          <word>;
        AMIMethod           <word>;
        coordinateSystem    <coordinateSystem>;
        fsm                 <bool>;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                            | Type  | Reqd | Deflt
      type     | Type name: turbulentDigitalFilterInlet | word  | yes  | -
      n        | Number of cells on turbulence generation plane       <!--
               --> (width height) or (e2 e3) | Vector2D\<label\> | yes | -
      mean     | Mean quantity             | PatchFunction1\<vector\> <!--
               -->                  or PatchFunction1\<scalar\> | yes  | -
      R        | Reynolds stresses     | PatchFunction1\<symmTensor\> <!--
               -->                  or PatchFunction1\<scalar\> | yes  | -
      L        | Integral scales            | tensor or vector  | yes  | -
      kernelType | Autocorrelation function form  | word | no   | Gaussian
      AMIMethod  | Mapping method            | word | no | faceAreaWeightAMI
      coordinateSystem | Coordinate system of turbulence generation plane <!--
               -->                           | coordinateSystem | no   | -
      fsm      | Flag to turn on the forward-stepwise method | bool | no | false
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink
      - \link PatchFunction1.H \endlink
      - \link AMIPatchToPatchInterpolation.H \endlink
      - \link coordinateSystem.H \endlink
      - \link IntegralScaleBox.H \endlink

    Options for the \c kernelType entry:
    \verbatim
      Gaussian    | Standard Gaussian function
      exponential | Exponential function
    \endverbatim

Note
  - The order of elements of the Reynolds stresses for the vector-based
    condition is "(Rxx Rxy Rxz Ryy Ryz Rzz)" and for the scalar-based
    condition is "R".
  - The order of elements of the integral scales for the vector-based
    condition is "(Lxu Lxv Lxw Lyu Lyv Lyw Lzu Lzv Lzw)" and for the
    scalar-based condition is "(Lxi Lyi Lzi)", where \c x here indicates
    the streamwise components and \c y/z lateral components. The streamwise
    components (i.e. \c x) are Lagrangian time scales in units of seconds
    and the remaining components are length scales in units of meters.
  - \c adjustTimeStep=true option is only supported by
    the forward-stepwise method.
  - The default global Cartesian coordinate system "(e1 e2 e3)"
    corresponds to "(x y z)".

See also
  - turbulentDFSEMInletFvPatchVectorField.C

SourceFiles
    turbulentDigitalFilterInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_turbulentDigitalFilterInletFvPatchField_H
#define Foam_turbulentDigitalFilterInletFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "PatchFunction1.H"
#include "IntegralScaleBox.H"
#include "AMIPatchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class turbulentDigitalFilterInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class turbulentDigitalFilterInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Typedefs

        //- Type of input set of Reynolds stresses
        typedef typename std::conditional
        <
            std::is_same<Type, scalar>::value,
            scalar,
            symmTensor
        >::type TypeR;


    // Private Data

        //- Pointer to AMI interpolator
        mutable autoPtr<AMIPatchToPatchInterpolation> AMIPtr_;

        //- Mean field on patch in global coordinates
        autoPtr<PatchFunction1<Type>> meanPtr_;

        //- Reynolds stresses field on patch in global coordinates
        autoPtr<PatchFunction1<TypeR>> Rptr_;

        //- Current time index
        label curTimeIndex_;

        //- Patch normal into the domain
        vector patchNormal_;

        //- Integral-scale field container in global coordinates
        turbulence::IntegralScaleBox<Type> L_;


    // Private Member Functions

        //- Return the patch normal
        vector calcPatchNormal() const;

        //- Initialise the patch content
        void initialisePatch();

        //- Map two-point correlations (integral scales)
        void mapL(Field<Type>& fld);

        //- Map one-point correlations (Reynolds stresses) for scalar fields
        void mapR(scalarField& fld) const;

        //- Map one-point correlations (Reynolds stresses) for vector fields
        void mapR(vectorField& fld) const;

        //- Map mean field for scalar fields
        void mapMean(scalarField& fld) const;

        //- Map mean field for vector fields
        void mapMean(vectorField& fld) const;


public:

   //- Runtime type information
   TypeName("turbulentDigitalFilterInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentDigitalFilterInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and patch field
        turbulentDigitalFilterInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDigitalFilterInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- turbulentDigitalFilterInletFvPatchField onto a new patch
        turbulentDigitalFilterInletFvPatchField
        (
            const turbulentDigitalFilterInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDigitalFilterInletFvPatchField
        (
            const turbulentDigitalFilterInletFvPatchField<Type>&
        );

        //- Construct as copy setting internal field reference
        turbulentDigitalFilterInletFvPatchField
        (
            const turbulentDigitalFilterInletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return fvPatchField<Type>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return fvPatchField<Type>::Clone(*this, iF);
        }


    // Member Functions

    // Mapping

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchField<Type>& ptf,
            const labelList& addr
        );


    // Evaluation

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


    // I-O

        //- Write boundary condition settings
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentDigitalFilterInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

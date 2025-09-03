/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::functionObjects::blendingFactor

Group
    grpFieldFunctionObjects

Description
    Computes blending factor as an indicator about which
    of the schemes is active across the domain.

    Blending factor values mean:
    \verbatim
      0   =  scheme 0
      1   =  scheme 1
      0-1 =  a blending factor between scheme 0 and scheme 1
    \endverbatim

    Operands:
    \table
      Operand       | Type           | Location
      input         | -              | -
      output file   | dat   | $FOAM_CASE/postProcessing/\<FO\>/\<time\>/\<file\>
      output field  | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    blendingFactor1
    {
        // Mandatory entries
        type            blendingFactor;
        libs            (fieldFunctionObjects);
        field           <field>;

        // Optional entries
        phi             phi;
        tolerance       0.001;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                           | Type   | Reqd | Deflt
      type      | Type name: blendingFactor             | word   | yes  | -
      libs      | Library name: fieldFunctionObjects    | word   | yes  | -
      field     | Name of the operand field             | word   | yes  | -
      phi       | Name of flux field                    | word   | no   | phi
      tolerance | Tolerance for number of blended cells | scalar | no   | 0.001
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink
     - \link writeFile.H \endlink

    Usage by the \c postProcess utility is not available.

Note
  - The \c blendingFactor function object overwrites the \c DEShybrid:Factor
    field internally when \c blendedSchemeBase debug flag is active.
    However, users are allowed to write out the original \c DEShybrid:Factor
    field by executing the \c writeObjects function object before
    any \c blendingFactor function object execution.

SourceFiles
    blendingFactor.C
    blendingFactorTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_blendingFactor_H
#define functionObjects_blendingFactor_H

#include "fieldExpression.H"
#include "writeFile.H"
#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class blendingFactor Declaration
\*---------------------------------------------------------------------------*/

class blendingFactor
:
    public fieldExpression,
    public writeFile
{
    // Private Data

        //- Name of flux field
        word phiName_;

        //- Tolerance used when calculating the number of blended cells
        scalar tolerance_;


    // Private Member Functions

        //- Calculate the blending factor field
        template<class Type>
        void calcBlendingFactor
        (
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const typename fv::convectionScheme<Type>& cs
        );

        //- Calculate the blending factor field
        template<class Type>
        bool calcScheme();

        //- Calculate the blending factor field and return true if successful
        virtual bool calc();


protected:

    // Protected Member Functions

        //- Write the file header
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("blendingFactor");


    // Constructors

        //- Construct from Time and dictionary
        blendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        blendingFactor(const blendingFactor&) = delete;

        //- No copy assignment
        void operator=(const blendingFactor&) = delete;


    //- Destructor
    virtual ~blendingFactor() = default;


    // Member Functions

        //- Read the blendingFactor data
        virtual bool read(const dictionary&);

        //- Write the blendingFactor
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blendingFactorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

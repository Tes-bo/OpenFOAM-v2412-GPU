/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Description
    Template for use with dynamic code generation of a
    fixedValue pointPatchField.

    - without state

SourceFiles
    fixedValuePointPatchFieldTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef coded_fixedValuePointPatchTemplate${FieldType}_H
#define coded_fixedValuePointPatchTemplate${FieldType}_H

#include "fixedValuePointPatchFields.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                A templated FixedValuePointPatchField
\*---------------------------------------------------------------------------*/

class ${typeName}FixedValuePointPatch${FieldType}
:
    public fixedValuePointPatchField<${TemplateType}>,
    public dictionaryContent
{
    //- The parent boundary condition type
    typedef fixedValuePointPatchField<${TemplateType}> parent_bctype;


    // Private Member Functions

        //- Report a message with the SHA1sum
        inline static void printMessage(const char* message)
        {
            Info<< message << " sha1: " << SHA1sum << '\n';
        }

public:

    //- SHA1 representation of the code content
    static constexpr const char* const SHA1sum = "${SHA1sum}";

    //- Runtime type information
    TypeName("${typeName}");


    // Constructors

        //- Construct from patch and internal field
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const dictionary& dict,
            IOobjectOption::readOption requireValue = IOobjectOption::MUST_READ
        );

        //- Construct by mapping a copy onto a new patch
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&,
            const pointPatch&,
            const DimensionedField<${TemplateType}, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&
        );

        //- Construct as copy setting internal field reference
        ${typeName}FixedValuePointPatch${FieldType}
        (
            const ${typeName}FixedValuePointPatch${FieldType}&,
            const DimensionedField<${TemplateType}, pointMesh>&
        );

        //- Return a clone
        virtual autoPtr<pointPatchField<${TemplateType}>> clone() const
        {
            return pointPatchField<${TemplateType}>::Clone(*this);
        }

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<${TemplateType}>> clone
        (
            const DimensionedField<${TemplateType}, pointMesh>& iF
        ) const
        {
            return pointPatchField<${TemplateType}>::Clone(*this, iF);
        }


    //- Destructor
    virtual ~${typeName}FixedValuePointPatch${FieldType}();


    // Member Functions

        //- Code context as a dictionary
        const dictionary& codeContext() const noexcept
        {
            return dictionaryContent::dict();
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

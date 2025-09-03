/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2024 OpenCFD Ltd.
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
    Template for use with dynamic code generation of a PatchFunction1

SourceFiles
    codedPatchFunction1Template.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCode_codedPatchFunction1_${typeName}_${FieldType}_H
#define dynamicCode_codedPatchFunction1_${typeName}_${FieldType}_H

#include "PatchFunction1.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                        A coded version PatchFunction1
\*---------------------------------------------------------------------------*/

class ${typeName}PatchFunction1${FieldType}
:
    public PatchFunction1<${TemplateType}>,
    public dictionaryContent
{
    //- The parent PatchFunction1 type
    typedef PatchFunction1<${TemplateType}> parent_bctype;

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
        ${typeName}PatchFunction1${FieldType}
        (
            const polyPatch& pp,
            const word& redirectType,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct, resetting patch
        ${typeName}PatchFunction1${FieldType}
        (
            const ${typeName}PatchFunction1${FieldType}& rhs,
            const polyPatch& pp
        );

        //- Copy construct
        ${typeName}PatchFunction1${FieldType}
        (
            const ${typeName}PatchFunction1${FieldType}& rhs
        ) = default;

        //- Return a clone
        virtual tmp<PatchFunction1<${TemplateType}>> clone() const
        {
            return PatchFunction1<${TemplateType}>::Clone(*this);
        }

        //- Return a clone, setting the patch
        virtual tmp<PatchFunction1<${TemplateType}>> clone
        (
            const polyPatch& pp
        ) const
        {
            return PatchFunction1<${TemplateType}>::Clone(*this, pp);
        }


    //- Destructor
    virtual ~${typeName}PatchFunction1${FieldType}() = default;


    // Member Functions

        //- Is value uniform (i.e. independent of coordinate)
        virtual bool uniform() const { return false; }

        //- Return value as a function of (scalar) independent variable
        virtual tmp<Field<${TemplateType}>> value(const scalar x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

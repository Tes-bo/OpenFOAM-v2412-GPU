/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Template for use with dynamic code generation of a Function1

SourceFiles
    codedFunction1Template.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCode_codedFunction1_${typeName}_${TemplateType}_H
#define dynamicCode_codedFunction1_${typeName}_${TemplateType}_H

#include "Function1.H"
#include "dictionaryContent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                           A coded version Function1
\*---------------------------------------------------------------------------*/

class ${typeName}Function1_${TemplateType}
:
    public Function1<${TemplateType}>,
    public dictionaryContent
{
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

        //- Construct from entry name, dictionary and registry
        ${typeName}Function1_${TemplateType}
        (
            const word& entryName,
            const dictionary& dict,
            const objectRegistry* obrPtr = nullptr
        );

        //- Copy construct
        ${typeName}Function1_${TemplateType}
        (
            const ${typeName}Function1_${TemplateType}& rhs
        ) = default;

        //- Return a clone
        virtual tmp<Function1<${TemplateType}>> clone() const
        {
            return Function1<${TemplateType}>::Clone(*this);
        }


    //- Destructor
    virtual ~${typeName}Function1_${TemplateType}() = default;


    // Member Functions

        //- Return value as a function of (scalar) independent variable
        virtual ${TemplateType} value(const scalar x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

Class
    Foam::fv::codedSource

Group
    grpFvOptionsSources

Description
    Constructs on-the-fly fvOption source

    The hook functions take the following arguments:

    \verbatim
    codeCorrect
    (
        GeometricField<Type, fvPatchField, volMesh>& field
    )

    codeAddSup
    (
        fvMatrix<Type>& eqn,
        const label fieldi
    )

    codeAddSupRho
    (
        const volScalarField& rho,
        fvMatrix<Type>& eqn,
        const label fieldi
    )

    codeConstrain
    (
        fvMatrix<Type>& eqn,
        const label fieldi
    )
    \endverbatim

    where :
        field is the name of the field in the fields list
        eqn is the fvMatrix

    These are in addition to the usual code entries:
    \plaintable
       codeInclude | include files
       codeOptions | compiler line: added to EXE_INC (Make/options)
       codeLibs    | linker line: added to LIB_LIBS (Make/options)
       localCode   | c++; local static functions
    \endplaintable

Usage
    Example usage in controlDict:
    \verbatim
    energySource
    {
        type            scalarCodedSource;

        scalarCodedSourceCoeffs
        {
            selectionMode   all;

            fields          (h);
            name            sourceTime;

            codeInclude
            #{
            #};

            codeCorrect
            #{
                Pout<< "**codeCorrect**" << endl;
            #};

            codeAddSup
            #{
                const Time& time = mesh().time();
                const scalarField& V = mesh_.V();
                scalarField& heSource = eqn.source();
                heSource -= 0.1*sqr(time.value())*V;
            #};

            codeConstrain
            #{
                Pout<< "**codeConstrain**" << endl;
            #};
        }
    }
    \endverbatim

SourceFiles
    codedSource.C

\*---------------------------------------------------------------------------*/

#ifndef fv_CodedFvSource_H
#define fv_CodedFvSource_H

#include "cellSetOption.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class codedSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedSource
:
    public fv::cellSetOption,
    protected codedBase
{
protected:

    // Protected Data

        word name_;

        string codeCorrect_;
        string codeAddSup_;
        string codeAddSupRho_;
        string codeConstrain_;

        //- Underlying code
        mutable autoPtr<fv::option> redirectOptionPtr_;


    // Protected Member Functions

        //- Mutable access to the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Description (type + name) for the output
        virtual string description() const;

        //- Clear redirected object(s)
        virtual void clearRedirect() const;

        //- Get the dictionary to initialize the code context
        virtual const dictionary& codeDict() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "codedFvOptionTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "codedFvOptionTemplate.H";


    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from components
        CodedSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Dynamically compiled fvOption
        fv::option& redirectOption() const;


    // Evaluation

        //- Correct field
        virtual void correct
        (
            GeometricField<Type, fvPatchField, volMesh>&
        );

        //- Explicit/implicit matrix contributions
        virtual void addSup
        (
            fvMatrix<Type>& eqn,
            const label fieldi
        );

        //- Explicit/implicit matrix contributions to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<Type>& eqn,
            const label fieldi
        );

        //- Set value
        virtual void constrain
        (
            fvMatrix<Type>& eqn,
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

#ifdef NoRepository
    #include "CodedFvSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

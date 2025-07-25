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
    Foam::fa::externalFileSource

Group
    grpFaOptionsSources

Description
    Applies sources on a specified field within a specified region
    by using an external table file.

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    externalFileSource1
    {
        // Mandatory entries (unmodifiable)
        type                externalFileSource;
        fieldName           <fieldName>;
        tableName           <tableFileName.dat>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type  | Reqd | Dflt
      type      | Type name: externalFileSource  | word  | yes  | -
      fieldName | Name of operand field          | word  | yes  | -
      tableName | Name of operand table file     | word  | yes  | -
      store     | Register external field 'pExt' | bool  | no   | false
    \endtable

    The inherited entries are elaborated in:
     - \link faOption.H \endlink
     - \link faceSetOption.H \endlink

See also
  - Foam::PatchFunction1Types

SourceFiles
    externalFileSource.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_externalFileSource_H
#define Foam_fa_externalFileSource_H

#include "faOption.H"
#include "areaFields.H"
#include "faceSetOption.H"
#include "MappedFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                        Class externalFileSource Declaration
\*---------------------------------------------------------------------------*/

class externalFileSource
:
    public fa::faceSetOption
{
    // Private Data

        //- Name of the field to apply this source
        word fieldName_;

        //- Name of the table
        word tableName_;

        //- External pressure field
        DimensionedField<scalar, areaMesh> pExt_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Mapped data from file
        PtrList<PatchFunction1Types::MappedFile<scalar>> mapping_;


    // Private Member Functions

        //- Update the pExt_ mapping
        void updateMapping();


public:

    //- Runtime type information
    TypeName("externalFileSource");


    // Constructors

        //- Construct from explicit source name and mesh
        externalFileSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        externalFileSource(const externalFileSource&) = delete;

        //- No copy assignment
        void operator=(const externalFileSource&) = delete;


    //- Destructor
    virtual ~externalFileSource() = default;


    // Member Functions

        // Evaluation

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const areaScalarField& rho,
                faMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

Namespace
    Foam::histogramModels

Description
    A namespace for various histogram model implementations.

Class
    Foam::histogramModel

Description
    A base class for histogram models.

SourceFiles
    histogramModel.C
    histogramModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_histogramModel_H
#define Foam_histogramModel_H

#include "writeFile.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class histogramModel Declaration
\*---------------------------------------------------------------------------*/

class histogramModel
:
    public functionObjects::writeFile
{
    // Private Data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Name of operand field
        word fieldName_;


protected:

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);

        //- Return requested field from the object registry
        //- or read+register the field to the object registry
        volScalarField& getOrReadField(const word& fieldName) const;

        //- Write histogram data
        void write
        (
            scalarField& dataNormalised,
            const labelField& dataCount,
            const scalarField& magMidBin
        );


public:

    //- Runtime type information
    TypeName("histogramModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            histogramModel,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected histogram model
        static autoPtr<histogramModel> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        histogramModel
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        histogramModel(const histogramModel&) = delete;

        //- No copy assignment
        void operator=(const histogramModel&) = delete;


    //- Destructor
    virtual ~histogramModel() = default;


    // Member Functions

    // Access

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return const reference to the operand field name
        const word& fieldName() const noexcept
        {
            return fieldName_;
        }


    // I-O

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Write data to stream and files
        virtual bool write(const bool log) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

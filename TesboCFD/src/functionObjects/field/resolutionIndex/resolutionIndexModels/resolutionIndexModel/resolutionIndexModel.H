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
    Foam::resolutionIndexModels

Description
    A namespace for various \c resolutionIndex model implementations.

Class
    Foam::resolutionIndexModel

Description
    A base class for \c resolutionIndex models.

SourceFiles
    resolutionIndexModel.C
    resolutionIndexModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_resolutionIndexModel_H
#define Foam_resolutionIndexModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class resolutionIndexModel Declaration
\*---------------------------------------------------------------------------*/

class resolutionIndexModel
{
    // Private Data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Name of result field
        word resultName_;


protected:

    // Protected Member Functions

        //- Return requested field from the object registry
        //- or read+register the field to the object registry
        template<class GeoFieldType>
        GeoFieldType& getOrReadField(const word& fieldName) const;

        //- Return cell volume field
        tmp<volScalarField> V() const;


public:

    //- Runtime type information
    TypeName("resolutionIndexModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            resolutionIndexModel,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Selectors

        //- Return a reference to the selected resolutionIndex model
        static autoPtr<resolutionIndexModel> New
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Constructors

        //- Construct from components
        resolutionIndexModel
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- No copy construct
        resolutionIndexModel(const resolutionIndexModel&) = delete;

        //- No copy assignment
        void operator=(const resolutionIndexModel&) = delete;


    //- Destructor
    virtual ~resolutionIndexModel() = default;


    // Member Functions

    // Access

        //- Return const reference to the mesh
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return const reference to the result name
        const word& resultName() const noexcept
        {
            return resultName_;
        }


    // I-O

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate the result field
        virtual bool execute() = 0;

        //- Write the result field
        virtual bool write() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "resolutionIndexModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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
    Foam::heatExchangerModels

Description
    A namespace for various heat exchanger model implementations.

Class
    Foam::heatExchangerModel

Description
    Base class for heat exchanger models to handle various
    characteristics for the \c heatExchangerSource fvOption.

SourceFiles
    heatExchangerModel.C
    heatExchangerModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_heatExchangerModel_H
#define Foam_heatExchangerModel_H

#include "fvMesh.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class heatExchangerModel Declaration
\*---------------------------------------------------------------------------*/

class heatExchangerModel
:
    public functionObjects::writeFile
{
protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionary containing coefficients specific to the chosen model
        const dictionary& coeffs_;

        //- Reference to the name of the fvOption source
        const word& name_;

        //- Name of operand velocity field
        word UName_;

        //- Name of operand temperature field
        word TName_;

        //- Name of operand flux field
        word phiName_;

        //- Name of the faceZone at the heat exchanger inlet
        word faceZoneName_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch IDs per face
        labelList facePatchId_;

        //- List of +1/-1 representing face flip map (1 use as is, -1 negate)
        labelList faceSign_;


public:

    //- Runtime type information
    TypeName("heatExchangerModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatExchangerModel,
            dictionary,
            (
                const fvMesh& mesh,
                const word& name,
                const dictionary& coeffs
            ),
            (mesh, name, coeffs)
        );


    // Selectors

        //- Return a reference to the selected heat exchanger model
        static autoPtr<heatExchangerModel> New
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& coeffs
        );


    // Constructors

        //- Construct from components
        heatExchangerModel
        (
            const fvMesh& mesh,
            const word& name,
            const dictionary& coeffs
        );

        //- No copy construct
        heatExchangerModel(const heatExchangerModel&) = delete;

        //- No copy assignment
        void operator=(const heatExchangerModel&) = delete;


    //- Destructor
    virtual ~heatExchangerModel() = default;


    // Member Functions

    // Access

        //- Return const reference to the name of velocity field
        virtual const word& U() const
        {
            return UName_;
        }


    // Evaluation

        //- Initialise data members of the model
        virtual void initialise();

        //- Return energy density per unit length [J/m3/m]
        virtual tmp<scalarField> energyDensity(const labelList& cells) = 0;


    // I-O

        //- Read top-level dictionary
        virtual bool read(const dictionary& dict) = 0;

        //- Write data to stream and files
        virtual void write(const bool log) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::VF

Description
    A namespace for various \c viewFactor model implementations.

Class
    Foam::VF::viewFactorModel

Description
    A base class for \c viewFactor models.

Usage
    Minimal example in \c <constant>/viewFactorsDict:
    \verbatim
    // Mandatory entries
    writeViewFactors    <bool>;

    // Optional entries
    writeRays           <bool>;
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                       | Type | Reqd | Deflt
      writeViewFactors | Flag to write the view factor field | bool | yes  | -
      writeRays   | Flag to write the ray geometry    | bool | no   | false
    \endtable

SourceFiles
    viewFactorModel.C
    viewFactorModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vf_viewFactorModel_H
#define Foam_vf_viewFactorModel_H

#include "autoPtr.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

namespace VF
{

// Forward Declarations
class raySearchEngine;

/*---------------------------------------------------------------------------*\
                       Class viewFactorModel Declaration
\*---------------------------------------------------------------------------*/

class viewFactorModel
{
protected:

    // Protected Data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Run-time selectable ray search engine
        autoPtr<raySearchEngine> searchEnginePtr_;

        //- Flag to write the view factor field
        bool writeViewFactors_;

        //- Flag to write the ray geometry
        bool writeRays_;


    // Protected Member Functions

        //- Write ray geometry to file
        static void writeRays
        (
            const fileName& fName,
            const pointField& compactCf,
            const labelListList& visibleFaceFaces
        );

        //- Calculate the view factors using run-time selectable model
        virtual scalarListList calculate
        (
            const labelListList& visibleFaceFaces,
            const pointField& compactCoarseCf,
            const vectorField& compactCoarseSf,
            const List<List<vector>>& compactFineSf,
            const List<List<point>>& compactFineCf,
            const DynamicList<List<point>>& compactPoints,
            const DynamicList<label>& compactPatchId
        ) const = 0;


public:

    //- Runtime type information
    TypeName("viewFactorModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        viewFactorModel,
        mesh,
        (
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (mesh, dict)
    );

    //- Selector
    static autoPtr<viewFactorModel> New
    (
        const fvMesh& mesh,
        const dictionary& dict
    );


    // Generated Methods

        //- No copy construct
        viewFactorModel(const viewFactorModel&) = delete;

        //- No copy assignment
        void operator=(const viewFactorModel&) = delete;


    //- Constructor
    viewFactorModel(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~viewFactorModel();


    // Public Member Functions

        //- Calculate the view factors
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace VF
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

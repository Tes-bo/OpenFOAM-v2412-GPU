/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::porosityModel

Description
    Top level model for porosity models

SourceFiles
    porosityModel.C
    porosityModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_porosityModel_H
#define Foam_porosityModel_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class porosityModel Declaration
\*---------------------------------------------------------------------------*/

class porosityModel
:
    public regIOobject
{
protected:

    // Protected Data

        //- Porosity name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dictionary used for model construction
        const dictionary dict_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Porosity active flag
        bool active_;

        //- Name(s) of cell-zone
        wordRe zoneName_;

        //- Cell zone IDs
        labelList cellZoneIDs_;

        //- Local coordinate system
        autoPtr<coordinateSystem> csysPtr_;


    // Protected Member Functions

        //- Transform the model data wrt mesh changes
        virtual void calcTransformModelData() = 0;

        //- Adjust negative resistance values to be multiplier of max value
        void adjustNegativeResistance(dimensionedVector& resist);

        //- Calculate the porosity force
        virtual void calcForce
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu,
            vectorField& force
        ) const = 0;

        virtual void correct(fvVectorMatrix& UEqn) const = 0;

        virtual void correct
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        ) const = 0;

        virtual void correct
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU
        ) const = 0;


        //- Local coordinate system
        inline const coordinateSystem& csys() const;

        //- Return label index
        inline label fieldIndex(const label index) const;

        //- No copy construct
        porosityModel(const porosityModel&) = delete;

        //- No copy assignment
        void operator=(const porosityModel&) = delete;


public:

    //- Runtime type information
    TypeName("porosityModel");

    //- Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        porosityModel,
        mesh,
        (
            const word& modelName,
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const wordRe& cellZoneName
        ),
        (modelName, name, mesh, dict, cellZoneName)
    );

    //- Constructor
    porosityModel
    (
        const word& name,
        const word& modelType,
        const fvMesh& mesh,
        const dictionary& dict,
        const wordRe& cellZoneName = wordRe::null
    );


    //- Return pointer to new porosityModel object created on the freestore
    //  from an Istream
    class iNew
    {
        //- Reference to the mesh database
        const fvMesh& mesh_;
        const word& name_;

    public:

        iNew
        (
            const fvMesh& mesh,
            const word& name
        )
        :
            mesh_(mesh),
            name_(name)
        {}

        autoPtr<porosityModel> operator()(Istream& is) const
        {
            const dictionary dict(is);

            return autoPtr<porosityModel>
            (
                porosityModel::New
                (
                    name_,
                    mesh_,
                    dict
                )
            );
        }
    };

    //- Selector
    static autoPtr<porosityModel> New
    (
        const word& name,
        const fvMesh& mesh,
        const dictionary& dict,
        const wordRe& cellZoneName = wordRe::null
    );


    //- Destructor
    virtual ~porosityModel() = default;


    // Member Functions

        //- Return const access to the porosity model name
        inline const word& name() const;

        //- Return const access to the porosity active flag
        inline bool active() const;

        //- Return const access to the cell zone IDs
        inline const labelList& cellZoneIDs() const;

        //- Return dictionary used for model construction
        const dictionary& dict() const;

        //- Transform the model data wrt mesh changes
        virtual void transformModelData();

        //- Return the force over the cell zone(s)
        virtual tmp<vectorField> force
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        virtual void addResistance(fvVectorMatrix& UEqn);

        //- Add resistance
        virtual void addResistance
        (
            fvVectorMatrix& UEqn,
            const volScalarField& rho,
            const volScalarField& mu
        );

        //- Add resistance
        virtual void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC
        );


    // I-O

        //- Write
        virtual bool writeData(Ostream& os) const;

        //- Inherit read from regIOobject
        using regIOobject::read;

        //- Is object global
        virtual bool global() const
        {
            return true;
        }

        //- Return complete path + object name if the file exists
        //  either in the case/processor or case otherwise null
        virtual fileName filePath() const
        {
            return globalFilePath(type());
        }

        //- Read porosity dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global file type for porosityModel
template<>
struct is_globalIOobject<porosityModel> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "porosityModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

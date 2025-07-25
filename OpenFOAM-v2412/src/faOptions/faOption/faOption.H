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
    Foam::fa::option

Description
    Base abstract class for handling finite area options (i.e. \c faOption).

Usage
    Minimal example by using \c constant/faOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory entries (unmodifiable)
        type              <faOptionName>;

        // Mandatory entries (runtime modifiable)
        region            <regionName>;

        // Optional entries (unmodifiable/runtime modifiable)
        <faOption>Coeffs
        {
            // subdictionary entries
        }

        // Optional entries (runtime modifiable)
        active            true;
        log               true;
    }
    \endverbatim

    where the entries mean:
    \table
      Property         | Description                   | Type  | Reqd | Dflt
      type             | Name of operand faOption      | word  | yes  | -
      region           | Name of operand region        | word  | yes  | -
      \<faOption\>Coeffs | Dictionary containing settings of <!--
                    --> the selected faOption settings | dictionary | no | -
      active   | Flag to (de)activate faOption         | bool  | no   | true
      log      | Flag to log faOption-related info     | bool  | no   | true
    \endtable

Note
  - Source/sink options are to be added to the right-hand side of equations.

SourceFiles
    faOption.C
    faOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fa_option_H
#define Foam_fa_option_H

#include "faMatricesFwd.H"
#include "areaFieldsFwd.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "volSurfaceMapping.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class option
{
    // Private Member Functions

        //- Construct region mesh and fields
        void constructMeshObjects();


protected:

    // Protected Data

        //- Source name
        const word name_;

        //- Model type
        const word modelType_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Top level source dictionary
        dictionary dict_;

        //- Dictionary containing source coefficients
        dictionary coeffs_;

        //- Field names to apply source to - populated by derived models
        wordList fieldNames_;

        //- Applied flag list - corresponds to each fieldNames_ entry
        List<bool> applied_;

        //- Region name (finite-area)
        word regionName_;


    // Protected Member Functions

        //- Resize/reset applied flag list for all fieldNames_ entries
        void resetApplied();


private:

    // Private Data

        //- Demand-driven: pointer to region mesh database
        mutable autoPtr<faMesh> regionMeshPtr_;

        //- Demand-driven: volume-to-surface mapping
        mutable autoPtr<volSurfaceMapping> vsmPtr_;

        //- Source active flag
        bool active_;


public:

    //- Switch write log to Info
    bool log;


    //- Runtime type information
    TypeName("option");


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            option,
            dictionary,
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, modelType, dict, mesh)
        );


    // Constructors

        //- Construct from components
        option
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<option> clone() const
        {
            NotImplemented;
            return nullptr;
        }

        //- Return pointer to new faOption object created
        //- on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh
            const fvMesh& mesh_;

            //- Name
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

            autoPtr<option> operator()(Istream& is) const
            {
                const dictionary dict(is);

                return autoPtr<option>
                (
                    option::New(name_, dict, mesh_)
                );
            }
        };


    // Selectors

        //- Return a reference to the selected faOption model
        static autoPtr<option> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~option() = default;


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const noexcept;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const noexcept;

            //- Return dictionary
            inline const dictionary& coeffs() const noexcept;

            //- True if source is active
            inline bool active() const noexcept;

            //- Set the applied flag to true for field index fieldi
            inline void setApplied(const label fieldi);

            //- The region name
            inline const word& regionName() const noexcept;

            //- Return the region mesh database (demand-driven)
            inline const faMesh& regionMesh() const;

            //- Return volSurfaceMapping (demand-driven)
            inline const volSurfaceMapping& vsm() const;


        // Edit

            //- Change source active flag, return previous value
            inline bool active(const bool on) noexcept;


        // Checks

            //- Is the source active?
            virtual bool isActive();

            //- Return index of field name if found in fieldNames list
            virtual label applyToField(const word& fieldName) const;

            //- Check that the source has been applied
            virtual void checkApplied() const;


        // Evaluation

            // Explicit and implicit sources

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    faMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Explicit and implicit sources for compressible equations

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const areaScalarField& h,
                    const areaScalarField& rho,
                    faMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Constraints

                virtual void constrain
                (
                    faMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    faMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Correction

                virtual void correct(areaScalarField& field);
                virtual void correct(areaVectorField& field);
                virtual void correct(areaSphericalTensorField& field);
                virtual void correct(areaSymmTensorField& field);
                virtual void correct(areaTensorField& field);


        // IO

            //- Write the source header information
            virtual void writeHeader(Ostream&) const;

            //- Write the source footer information
            virtual void writeFooter(Ostream&) const;

            //- Write the source properties
            virtual void writeData(Ostream&) const;

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
    Copyright (C) 2020,2023 PCOpt/NTUA
    Copyright (C) 2020,2023 FOSS GP
------------------------------------------------------------------------------
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
    Foam::fv::option

Description
    Base abstract class for handling finite volume options (i.e. \c fvOption).

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    <userDefinedName1>
    {
        // Mandatory entries (unmodifiable)
        type              <fvOptionName>;

        // Optional entries (unmodifiable/runtime modifiable)
        <fvOption>Coeffs
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
      type             | Name of operand fvOption      | word  | yes  | -
      \<fvOption\>Coeffs | Dictionary containing settings of <!--
                    --> the selected fvOption settings | dictionary | no | -
      active   | Flag to (de)activate fvOption         | bool  | no   | true
      log      | Flag to log fvOption-related info     | bool  | no   | true
    \endtable

SourceFiles
    fvOption.C
    fvOptionIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvOption_H
#define Foam_fvOption_H

#include "fvMatricesFwd.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class option Declaration
\*---------------------------------------------------------------------------*/

class option
{
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

        //- Source active flag
        bool active_;


    // Protected Member Functions

        //- Resize/reset applied flag list for all fieldNames_ entries
        void resetApplied();


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

        //- Return pointer to new fvOption object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh
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

        //- Return a reference to the selected fvOption model
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
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Explicit and implicit sources for compressible equations

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Explicit and implicit sources for phase equations

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void addSup
                (
                    const volScalarField& alpha,
                    const volScalarField& rho,
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Constraints

                virtual void constrain
                (
                    fvMatrix<scalar>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<vector>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<sphericalTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<symmTensor>& eqn,
                    const label fieldi
                );

                virtual void constrain
                (
                    fvMatrix<tensor>& eqn,
                    const label fieldi
                );


            // Correction

                virtual void correct(volScalarField& field);
                virtual void correct(volVectorField& field);
                virtual void correct(volSphericalTensorField& field);
                virtual void correct(volSymmTensorField& field);
                virtual void correct(volTensorField& field);

                virtual void correct(surfaceScalarField& field);
                virtual void correct(surfaceVectorField& field);
                virtual void correct(surfaceTensorField& field);


            // Post process sensitivity field related to the fvOption

                virtual void postProcessSens
                (
                    scalarField& sensField,
                    const word& fieldName = word::null,
                    const word& designVariablesName = word::null
                );
                virtual void postProcessSens
                (
                    vectorField& sensField,
                    const word& fieldName = word::null,
                    const word& designVariablesName = word::null
                );
                virtual void postProcessSens
                (
                    tensorField& sensField,
                    const word& fieldName = word::null,
                    const word& designVariablesName = word::null
                );

                virtual void postProcessAuxSens
                (
                    const volScalarField& primalField,
                    const volScalarField& adjointField,
                    scalarField& sensField,
                    const word& fieldName = word::null
                );
                virtual void postProcessAuxSens
                (
                    const volVectorField& primalField,
                    const volVectorField& adjointField,
                    scalarField& sensField,
                    const word& fieldName = word::null
                );
                virtual void postProcessAuxSens
                (
                    const volTensorField& primalField,
                    const volTensorField& adjointField,
                    scalarField& sensField,
                    const word& fieldName = word::null
                );


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

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvOptionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

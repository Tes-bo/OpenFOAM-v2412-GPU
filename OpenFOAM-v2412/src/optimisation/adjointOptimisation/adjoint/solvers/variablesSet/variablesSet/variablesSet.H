/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::variablesSet

Description
    Base class for creating a set of variables

\*---------------------------------------------------------------------------*/

#ifndef variablesSet_H
#define variablesSet_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class variablesSet Declaration
\*---------------------------------------------------------------------------*/

class variablesSet
{
private:

    // Private Member Functions

        //- No copy construct
        variablesSet(const variablesSet&) = delete;

        //- No copy assignment
        void operator=(const variablesSet&) = delete;

        //- Add the solverName to every patch entry read from the boundaryField
        //- of the IOobject and allocate the corresponding GeometricField.
        //- Used to pass the adjoint solver name to the adjoint boundary
        //- conditions.
        //  Returns naked pointer but used only to feed autoPtrs and tmps,
        //  so memory handling shouldn't be an issue
        template<class Type, template<class> class PatchField, class GeoMesh>
        static GeometricField<Type, PatchField, GeoMesh>* allocateNamedField
        (
            const fvMesh& mesh,
            const IOobject& io,
            const word& solverName
        );

        //- Read field with base or custom field name
        template<class Type, template<class> class PatchField, class GeoMesh>
        static bool readFieldOK
        (
            autoPtr<GeometricField<Type, PatchField, GeoMesh>>& fieldPtr,
            const fvMesh& mesh,
            const word& baseName,
            const word& solverName,
            const bool useSolverNameForFields
        );


protected:

    // Protected data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- Solver name owning the variables set
        word solverName_;

        //- Append the solver name to the variables names?
        bool useSolverNameForFields_;


    // Protected Member Functions

        template<class Type, template<class> class PatchField, class GeoMesh>
        autoPtr<GeometricField<Type, PatchField, GeoMesh>> allocateRenamedField
        (
            const autoPtr<GeometricField<Type, PatchField, GeoMesh>>& bf
        );

        //- Swap autoPtrs and rename managed fields
        template<class Type, template<class> class PatchField, class GeoMesh>
        void swapAndRename
        (
            autoPtr<GeometricField<Type, PatchField, GeoMesh>>& p1,
            autoPtr<GeometricField<Type, PatchField, GeoMesh>>& p2
        );


public:


    // Static Data Members

        //- Run-time type information
        TypeName("variablesSet");


    // Constructors

        //- Construct from mesh and solver name
        variablesSet
        (
            fvMesh& mesh,
            const dictionary& dict
        );

        //- Clone the variablesSet
        //  Not implemented in base. Should be overloaded in derived classes
        virtual autoPtr<variablesSet> clone() const;


    //- Destructor
    virtual ~variablesSet() = default;


    // Member Functions

        // Access

            //- Return solver name
            const word& solverName() const;

            //- Append solver name to fields?
            bool useSolverNameForFields() const;

        // Set functions. Static in order to be used by other classes as well

            //- Read vol fields
            template<class Type>
            static void setField
            (
                autoPtr<GeometricField<Type, fvPatchField, volMesh>>& fieldPtr,
                const fvMesh& mesh,
                const word& baseName,
                const word& solverName,
                const bool useSolverNameForFields
            );

            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh>> allocateField
            (
                const fvMesh& mesh,
                const word& baseName,
                const word& solverName,
                const bool useSolverNameForFields
            );

            //- Turbulence model always reads fields with the prescribed name
            //- If a custom name is supplied, check whether this field exists,
            //- copy it to the field known by the turbulence model
            //- and re-name the latter
            template<class Type>
            void renameTurbulenceField
            (
                GeometricField<Type, fvPatchField, volMesh>& baseField,
                const word& solverName
            );

            //- Set flux field
            static void setFluxField
            (
                autoPtr<surfaceScalarField>& fieldPtr,
                const fvMesh& mesh,
                const volVectorField& velocity,
                const word& baseName,
                const word& solverName,
                const bool useSolverNameForFields
            );

            tmp<surfaceScalarField> allocateFluxField
            (
                const fvMesh& mesh,
                const volVectorField& velocity,
                const word& baseName,
                const word& solverName,
                const bool useSolverNameForFields
            );

            //- Auto create variable for mesh movement
            static tmp<volVectorField> autoCreateMeshMovementField
            (
                const fvMesh& mesh,
                const word& name,
                const dimensionSet& dims
            );

            //- Nullify field and old times, if present
            template<class Type, template<class> class PatchField, class GeoMesh>
            static void nullifyField
            (
                GeometricField<Type, PatchField, GeoMesh>& fieldPtr
            );

            //- Transfer the fields of another variablesSet to this
            virtual void transfer(variablesSet& vars);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "variablesSetTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

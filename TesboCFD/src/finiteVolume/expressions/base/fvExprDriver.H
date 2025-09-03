/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2010-2018 Bernhard Gschaider
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
    Foam::expressions::fvExprDriver

Description
    Base driver for parsing value expressions associated with an fvMesh.

    Largely based on code and ideas from swak4foam

    Properties
    \table
        Property     | Description                          | Required | Default
        variables    | List of variables for expressions    | no  | ()
        delayedVariables | List of delayed variables        | no  | ()
        storedVariables | List of stored variables          | no  | ()
        globalScopes | Scopes for global variables          | no  | ()
        allowShadowing | Allow variables to shadow field names | no  | false
    \endtable

    Debug Properties
    \table
        Property      | Description                         | Required | Default
        debug.driver  | Debug level (int) for base driver   | no  |
        debug.scanner | Add debug for scanner               | no  | false
        debug.parser  | Add debug for parser                | no  | false
    \endtable

SourceFiles
    fvExprDriverI.H
    fvExprDriver.C
    fvExprDriverIO.C
    fvExprDriverNew.C
    fvExprDriverTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_fvExprDriver_H
#define expressions_fvExprDriver_H

#include "exprDriver.H"
#include "exprResultDelayed.H"
#include "exprResultStored.H"
#include "pointMesh.H"
#include "volFields.H"
#include "topoSetSource.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

// Forward Declarations
class fvExprDriverWriter;

/*---------------------------------------------------------------------------*\
                         Class fvExprDriver Declaration
\*---------------------------------------------------------------------------*/

class fvExprDriver
:
    public expressions::exprDriver
{
    // Static Data

        //- Pointer to the "default" mesh
        static const fvMesh *defaultMeshPtr_;

        //- Cache cellSets, faceSets instead of reading from disc each time
        static bool cacheSets_;


    // Private Data

    // Stored Data

        //- The scopes for global variables
        List<word> globalScopes_;

        //- The (delayed) variables table
        HashTable<exprResultDelayed> delayedVariables_;

        //- Stored expressions. Read from dictionary and updated as required
        List<exprResultStored> storedVariables_;

        //- Time index when handling special variables
        label specialVariablesIndex_;

        //- The name of the other mesh (if it is to be required)
        word otherMeshName_;

        //- Writing and restoring
        autoPtr<fvExprDriverWriter> writer_;


    // Private Member Functions

        //- Read the IOobject and return the headerClassName
        static word getHeaderClassName
        (
            const polyMesh& mesh,
            const word& name
        );

        //- Read the set IOobject and return its headerClassName
        static word getSetClassName
        (
            const polyMesh& mesh,
            const word& name
        );


        //- No copy assignment
        void operator=(const fvExprDriver&) = delete;


protected:

    // Static Member Functions

        //- Determine mesh or region mesh as specified in the dictionary
        //- with the keyword "region"
        static const fvMesh& regionMesh
        (
            const dictionary& dict,
            const fvMesh& mesh,
            bool readIfNecessary
        );

        //- Default boundary type is calculated
        template<class T>
        static inline word defaultBoundaryType(const T&)
        {
            return "calculated";
        }

        //- Default boundary type for volume fields is zeroGradient
        //- since they are essentially just internal fields.
        template<class Type>
        static inline word defaultBoundaryType
        (
            const VolumeField<Type>&
        )
        {
            return "zeroGradient";
        }


        //- Apply correctBoundaryConditions (volume fields only)
        template<class T>
        static inline void correctField(T&) {}

        template<class Type>
        static inline void correctField
        (
            VolumeField<Type>& fld
        )
        {
            fld.correctBoundaryConditions();
        }


    // Protected Member Functions

    // Mesh related

        //- The mesh we are attached to
        virtual const fvMesh& mesh() const = 0;


    // Variables

        //- Define scopes for global variables
        void setGlobalScopes(const wordUList& scopes)
        {
            globalScopes_ = scopes;
        }

        //- Non-const access to the named variable (sub-classes only)
        inline virtual exprResult& variable(const word& name);

        //- Test existence of a global variable
        template<class T>
        bool isGlobalVariable
        (
            const word& name,
            const bool wantPointData = false,
            const label expectedSize = -1
        ) const;

        //- Return the global variable if available or a null result
        const exprResult& lookupGlobal(const word& name) const;


    // Fields

        //- Test for the existence of a mesh field
        template<class Type>
        bool isField
        (
            const word& name,
            const bool wantPointData = false,
            const label expectSize = -1  //!< ignored
        ) const;


        //- Retrieve field from memory or disk
        template<class GeomField>
        inline tmp<GeomField> getOrReadField
        (
            const word& name,
            const bool mandatory = true,
            const bool getOldTime = false
        );

        //- Retrieve point field from memory or disk
        template<class GeomField>
        inline tmp<GeomField> getOrReadPointField
        (
            const word& name,
            const bool mandatory = true,
            const bool getOldTime = false
        );

        //- Retrieve field from memory or disk (implementation)
        template<class GeomField, class MeshRef>
        tmp<GeomField> getOrReadFieldImpl
        (
            const word& name,
            const MeshRef& meshRef,
            const bool mandatory = true,
            const bool getOldTime = false
        );

        //- Helper function for getOrReadField
        template<class GeomField, class MeshRef>
        inline tmp<GeomField> readAndRegister
        (
            const word& name,
            const MeshRef& meshRef
        );


    // Sets

        //- The origin of the topoSet
        enum SetOrigin { INVALID = 0, NEW, FILE, MEMORY, CACHE };

        //- Get topoSet
        template<class T>
        autoPtr<T> getTopoSet
        (
            const fvMesh& mesh,
            const word& setName,
            SetOrigin& origin
        ) const;

        //- Update topoSet
        template<class T>
        inline bool updateSet
        (
            autoPtr<T>& setPtr,
            const word& setName,
            SetOrigin origin
        ) const;


    // Updating

        //- Examine current variable values and update stored variables
        virtual void updateSpecialVariables(bool force=false);

        //- Do we need a data file to be written
        virtual bool hasDataToWrite() const;

        //- Prepare/update special variables and add to dictionary,
        //- normally via the reader/writer
        virtual void prepareData(dictionary& dict) const;

        //- Read data from dictionary, normally via the reader/writer
        virtual void getData(const dictionary& dict);


public:

    // Friends
    friend class fvExprDriverWriter;


    // Static Member Functions

        //- Get the default mesh, if one is defined
        static const fvMesh& defaultMesh();

        //- Set the default mesh (if not already set)
        static const fvMesh* resetDefaultMesh
        (
            const fvMesh& mesh,
            const bool force = false  //!< Force reset, even if already set
        );


    //- Runtime type information
    TypeName("fvExprDriver");


    // Run-time selection

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvExprDriver,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict,mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvExprDriver,
            idName,
            (
                const word& ident,
                const fvMesh& mesh
            ),
            (ident, mesh)
        );


    // Constructors

        //- Default construct, and default construct with search preferences
        explicit fvExprDriver
        (
            enum exprDriver::searchControls search =
                exprDriver::searchControls::DEFAULT_SEARCH,
            const dictionary& dict = dictionary::null
        );

        //- Copy construct with dictionary reference
        explicit fvExprDriver
        (
            const fvExprDriver& rhs,
            const dictionary& dict
        );

        //- Construct from a dictionary
        explicit fvExprDriver(const dictionary& dict);


        //- Return a reference to the selected value driver
        static autoPtr<fvExprDriver> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return a reference to the selected value driver
        static autoPtr<fvExprDriver> New
        (
            const dictionary& dict
        );

        //- Return a reference to the selected value driver
        static autoPtr<fvExprDriver> New
        (
            const word& type,
            const word& id,
            const fvMesh& mesh
        );

        //- Not generally clonable
        virtual autoPtr<fvExprDriver> clone() = delete;


    //- Destructor
    virtual ~fvExprDriver();


    // Public Member Functions

        //- The natural field size for the expression
        virtual label size() const = 0;

        //- The point field size for the expression
        virtual label pointSize() const = 0;


    // General Controls

        //- Status of cache-sets (static variable)
        bool cacheSets() const { return cacheSets_; }


    // Variables

        //- Clear temporary variables, reset from expression strings
        virtual void clearVariables();

        //- True if named variable exists
        inline virtual bool hasVariable(const word& name) const;

        //- Return const-access to the named variable
        inline virtual const exprResult& variable(const word& name) const;

        //- Test for existence of a local/global variable or a field
        template<class Type>
        inline bool isVariableOrField
        (
            const word& name,
            const bool wantPointData = false,
            const label expectSize = -1
        ) const;

        //- Retrieve local/global variable as a tmp field
        //
        //  \param name The name of the local/global field
        //  \param expectSize  The size check on the variable, -1 to ignore
        //  \param mandatory A missing variable is Fatal, or return
        //      an invalid tmp
        template<class Type>
        tmp<Field<Type>> getVariable
        (
            const word& name,
            const label expectSize,
            const bool mandatory = true
        ) const;


        //- Lookup the field class name (memory or read from disk)
        //
        //  Return empty if the name cannot be resolved.
        word getFieldClassName(const word& name) const;


    // Types

        //- Return cell/face/point set type or unknown
        topoSetSource::sourceType topoSetType(const word& name) const;

        //- Return cell/face/point zone type or unknown
        topoSetSource::sourceType topoZoneType(const word& name) const;

        //- Return cell/face/point zone/set type or unknown
        topoSetSource::sourceType topoSourceType(const word& name) const;

        //- Get the labels associated with the topo set
        refPtr<labelList> getTopoSetLabels
        (
            const word& name,
            enum topoSetSource::sourceType setType
        ) const;

        //- Test if name is a known cellZone
        bool isCellZone(const word& name) const;

        //- Test if name is a known faceZone
        bool isFaceZone(const word& name) const;

        //- Test if name is a known pointZone
        bool isPointZone(const word& name) const;

        //- Test if name is a known cellSet
        bool isCellSet(const word& name) const;

        //- Test if name is a known faceSet
        bool isFaceSet(const word& name) const;

        //- Test if name is a known pointSet
        bool isPointSet(const word& name) const;


    // Evaluation

        //- Evaluate the expression
        //- and save as the specified named variable
        virtual void evaluateVariable
        (
            const word& varName,
            const expressions::exprString& expr
        );

        //- Evaluate an expression on a remote
        //- and save as the specified named variable
        //
        //  The fully qualified form of the remote is given as follows
        //  \verbatim
        //      type'name/region
        //  \endverbatim
        //
        //  If not specified, the default type is "patch", which means the
        //  following are equivalent
        //  \verbatim
        //      patch'name/region
        //      name/region
        //  \endverbatim
        //
        //  If region is identical to the current region, it can be omitted:
        //  \verbatim
        //      patch'name
        //      name   (default is patch)
        //  \endverbatim
        virtual void evaluateVariableRemote
        (
            string remote,
            const word& varName,
            const expressions::exprString& expr
        );


    // Fields

        //- Test existence of a local/global variable
        template<class Type>
        inline bool isVariable
        (
            const word& name,
            bool wantPointData = false,
            label expectSize = -1
        ) const;

        //- Test if specified field can be found in memory or disk
        template<class Type>
        bool foundField(const word& name) const;

        //- Read the IOobject for fieldName and return its headerClassName
        //  Empty if the field could not be found.
        word getTypeOfField(const word& fieldName) const;


    // Handling remote data (future)

        // //- Access the other mesh name (future)
        // const word& otherMeshName() const { return otherMeshName_; }
        //
        // //- Access the other mesh name (future)
        // word& otherMeshName() { return otherMeshName_; }


    // Reading

        //- Read variables, tables etc.
        //  Also usable for objects not constructed from a dictionary.
        virtual bool readDict(const dictionary& dict);


    // Writing

        //- Write "variables", "storedVariables", "delayedVariables",
        //- "globalScopes" if they are present.
        Ostream& writeCommon(Ostream& os, bool debug=false) const;

        //- Create a writer for this object
        void createWriterAndRead(const word& name);

        //- Write data if appropriate
        //- Should enable exact restarts
        void tryWrite() const;


    // Plugins (future)

        /// //- Tests for a plugin-function
        /// virtual bool hasPlugin(const word& name) = 0;
        ///
        /// //- Return a new plugin-function
        /// virtual autoPtr<CommonPluginFunction>
        /// getPlugin(const word& name) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvExprDriverI.H"

#ifdef NoRepository
    #include "fvExprDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

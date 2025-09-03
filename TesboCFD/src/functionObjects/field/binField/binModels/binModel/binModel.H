/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::binModel

Description
    Base class for bin models to handle general bin characteristics.

SourceFiles
    binModel.C
    binModelNew.C
    binModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_binModel_H
#define Foam_binModel_H

#include "dictionary.H"
#include "HashSet.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"
#include "coordinateSystem.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                            Class binModel Declaration
\*---------------------------------------------------------------------------*/

class binModel
:
    public functionObjects::writeFile
{
protected:

    // Protected Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Decompose patch values into normal and tangential components
        bool decomposePatchValues_;

        //- Flag to accumulate bin data with
        //- increasing distance in binning direction
        bool cumulative_;

        //- Local coordinate system of bins
        autoPtr<coordinateSystem> coordSysPtr_;

        //- Total number of bins
        label nBin_;

        //- Indices of operand patches
        labelList patchIDs_;

        //- Names of operand fields
        wordList fieldNames_;

        //- Indices of operand cell zones
        labelList cellZoneIDs_;

        //- List of file pointers; 1 file per field
        PtrList<OFstream> filePtrs_;


    // Protected Member Functions

        //- Set the co-ordinate system from dictionary and axes names
        void setCoordinateSystem
        (
            const dictionary& dict,
            const word& e3Name = word::null,
            const word& e1Name = word::null
        );

        //- Helper function to decompose patch values
        //- into normal and tangential components
        template<class Type>
        bool decomposePatchValues
        (
            List<List<Type>>& data,
            const label bini,
            const Type& v,
            const vector& n
        ) const;

        //- Helper function to construct a string description for a given type
        template<class Type>
        string writeComponents(const word& stem) const;

        //- Write binned data to stream
        template<class Type>
        void writeBinnedData
        (
            List<List<Type>>& data,
            Ostream& os
        ) const;


public:

    //- Runtime type information
    TypeName("binModel");


    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            binModel,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word& outputPrefix
            ),
            (dict, mesh, outputPrefix)
        );


    // Selectors

        //- Return a reference to the selected bin model
        static autoPtr<binModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& outputPrefix
        );


    // Constructors

        //- Construct from components
        binModel
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& outputPrefix
        );

        //- No copy construct
        binModel(const binModel&) = delete;

        //- No copy assignment
        void operator=(const binModel&) = delete;


    //- Destructor
    virtual ~binModel() = default;


    // Member Functions

        //- Read the dictionary
        virtual bool read(const dictionary& dict);


    // Access

        //- Return the total number of bins
        label nBin() const noexcept
        {
            return nBin_;
        };


    // Evaluation

        //- Initialise bin properties
        virtual void initialise() = 0;

        //- Apply bins
        virtual void apply() = 0;

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
bool binModel::decomposePatchValues
(
    List<List<vector>>& data,
    const label bini,
    const vector& v,
    const vector& n
) const;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "binModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

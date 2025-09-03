/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::catalyst::timeQuery

Description
    Low-level handling for Catalyst time queries

SourceFiles
    catalystTools.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catalystTools_H
#define functionObjects_catalystTools_H

#include "wordList.H"
#include "DynamicList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkCPDataDescription;
class vtkCPProcessor;
class vtkMultiBlockDataSet;
template<class T> class vtkSmartPointer;

namespace Foam
{

// Forward Declarations
class Ostream;
class Time;

namespace catalyst
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Container holding VTK outputs per Catalyst channel name
typedef HashTable<vtkSmartPointer<vtkMultiBlockDataSet>> outputChannels;

//- Print information about data description
void printInfo(Ostream& os, vtkCPDataDescription* descrip);


/*---------------------------------------------------------------------------*\
                     Class catalyst::timeQuery Declaration
\*---------------------------------------------------------------------------*/

//- Simple structure for Catalyst time queries
struct timeQuery
{
    double timeValue;
    long   timeIndex;
    bool   forced;

    constexpr timeQuery
    (
        double val,
        long index,
        bool forcedOutput=false
    ) noexcept
    :
        timeValue(val),
        timeIndex(index),
        forced(forcedOutput)
    {}

    timeQuery(const Foam::Time& currTime);
};


/*---------------------------------------------------------------------------*\
                     Class catalyst::dataQuery Declaration
\*---------------------------------------------------------------------------*/

//- Simple structure for Catalyst data description queries
class dataQuery
:
    public timeQuery
{
    //- Catalyst channel names for query, in the order of addition
    DynamicList<word> names_;

    //- Fields per channel
    HashTable<wordHashSet> fields_;

    //- Indices (into names_) of active channels
    labelHashSet active_;

    //- The data description
    vtkCPDataDescription* descrip_;

    //- Empty dummy hash of field names
    static const wordHashSet emptyWordHash_;


    //- No copy construct
    dataQuery(const dataQuery&) = delete;

    //- No copy assignment
    void operator=(const dataQuery&) = delete;


public:

    // Constructors

        //- Construct query for given time. Add channels later
        explicit dataQuery(const Foam::Time& when);


    //- Destructor
    ~dataQuery();


    // Member Functions

        //- Pointer to coprocess data description
        inline vtkCPDataDescription* description() const;

        //- Clear the channel information
        inline void clear();

        //- Empty if there are no active channels
        inline bool empty() const noexcept;

        //- The number of active channels
        inline label size() const noexcept;

        //- True if given channel exists (in fields)
        inline bool contains(const word& channel) const;

        //- The active channel numbers
        inline List<label> active() const;

        //- The (active) channel names
        inline wordList names() const;

        //- Fields for the given channel
        inline const wordHashSet& fields(const word& channel) const;

        //- Clear the specified channel
        inline void clear(const word& channel);

        //- Set fields for the specified channel and make it active
        void set(const word& channel, const wordHashSet& fields);

        //- Query coprocess, transform defined input channels to
        //- known requests.
        //  On input the dataQuery contains the channel names and their
        //  associated published fields.
        //  On output, retains the requested channels only.
        label query(vtkCPProcessor* coproc);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace catalyst

// Ostream Operator
Ostream& operator<<(Ostream& os, const catalyst::timeQuery& when);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "catalystToolsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

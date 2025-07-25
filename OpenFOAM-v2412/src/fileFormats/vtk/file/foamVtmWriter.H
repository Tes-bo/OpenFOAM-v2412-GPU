/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::vtk::vtmWriter

Description
    Provides a means of accumulating file entries for generating
    a vtkMultiBlockDataSet (.vtm) file.

    For example, to generate the following content:
    \verbatim
    <?xml version='1.0'?>
    <VTKFile type='vtkMultiBlockDataSet' ...>
      <vtkMultiBlockDataSet>
        <DataSet name='internal' file='internal.vtu' />
        <Block name='boundary'>
          <DataSet name='inlet' file='boundary/inlet.vtp' />
          <DataSet name='outlet' file='boundary/outlet.vtp' />
        </Block>
      </vtkMultiBlockDataSet>
      <FieldData>
        <DataArray type='Float32' Name='TimeValue' ...>
           12.345
        </DataArray>
      </FieldData>
    </VTKFile>
    \endverbatim

    The following code would be used:
    \code
        vtm.clear();
        vtm.setTime(12.345);

        vtm.append("internal", "internal.vtu");

        vtm.beginBlock("boundary");
        vtm.append("boundary/inlet.vtp");
        vtm.append("boundary/outlet.vtp");

        vtm.write("outputName");
    \endcode

SourceFiles
    foamVtmWriter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtmWriter_H
#define Foam_vtk_vtmWriter_H

#include "foamVtkOutputOptions.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class vtk::vtmWriter Declaration
\*---------------------------------------------------------------------------*/

class vtmWriter
{
    //- Simple structure for containing entries
    struct vtmEntry
    {
        enum Type
        {
            NONE = 0,
            DATA = 'D',
            BEGIN_BLOCK = '{',
            END_BLOCK = '}'
        };

        //- The entry type
        int type_;

        //- The 'name' entry (to describe block or data)
        string name_;

        //- The 'file' entry (data only)
        fileName file_;

        // Constructors

            vtmEntry(const vtmEntry&) = default;
            vtmEntry(vtmEntry&&) = default;
            vtmEntry& operator=(const vtmEntry&) = default;
            vtmEntry& operator=(vtmEntry&&) = default;

            //- Default construct
            vtmEntry()
            :
                type_(NONE)
            {}

            //- Construct from components
            vtmEntry(int what, const string& name, const fileName& file)
            :
                type_(what), name_(name), file_(file)
            {}


        // Factory Methods

            static vtmEntry block(const string& name)
            {
                return vtmEntry(BEGIN_BLOCK, name, "");
            }

            static vtmEntry endblock()
            {
                return vtmEntry(END_BLOCK, "", "");
            }

            static vtmEntry entry(const fileName& file)
            {
                return vtmEntry(DATA, "", file);
            }

            static vtmEntry entry(const string& name, const fileName& file)
            {
                return vtmEntry(DATA, name, file);
            }


        // Member Functions

            //- Test the type
            bool isType(Type what) const noexcept
            {
               return type_ == what;
            }

            //- Reset to NONE
            void clear();

            //- True if the entry is good.
            bool good() const noexcept;

            //- Output valid entry as XML
            bool write(vtk::formatter& format) const;
    };


    // Private Member Data

        //- Auto-generate names from 'file' entry?
        bool autoName_;

        //- Has a TimeValue for FieldData?
        bool hasTime_;

        //- A vtm file entry: begin/end block, dataset
        DynamicList<vtmEntry> entries_;

        //- LIFO stack of current block names
        DynamicList<word> blocks_;

        //- TimeValue for FieldData
        scalar timeValue_;


    // Private Member Functions

        //- Remove NONE entries
        bool pruneEmpty();

        //- Remove empty blocks
        bool pruneEmptyBlocks();

        //- Collapse block if it has a single dataset and the names allow it
        bool collapseBlocks();


public:

    // Constructors

        //- Default construct, with autoName on
        vtmWriter();

        //- Construct with specified behaviour for autoName
        explicit vtmWriter(bool autoName);


    //- Destructor
    ~vtmWriter() = default;


    // Member Functions

        //- File extension (always "vtm")
        inline static word ext();

        //- If there are no data sets
        bool empty() const;

        //- The number of data sets
        label size() const;


    // Content Management

        //- Clear all entries and reset output
        void clear();

        //- Define "TimeValue" for FieldData (name as per Catalyst output)
        void setTime(scalar timeValue);

        //- Define "TimeValue" for FieldData (name as per Catalyst output)
        void setTime(const Time& t);


        //- Start a new block, optionally with a name
        //  \return block depth
        label beginBlock(const word& blockName = word::null);

        //- End the previous block, optionally with name checking
        //  \return block depth
        label endBlock(const word& blockName = word::null);


        //- Add a file. The name is either empty or created with autoName
        //  \return True if file is non-empty
        bool append(const fileName& file);

        //- Add a file with name
        //  \return True if file is non-empty
        bool append(const word& name, const fileName& file);

        //- Add a file with given contentType extension
        //- The name is either empty or created with autoName
        //  \return True if file is non-empty
        bool append(const fileName& file, vtk::fileTag contentType);

        //- Add a file with name, with given contentType extension
        //  \return True if file is non-empty
        bool append
        (
            const word& name,
            const fileName& file,
            vtk::fileTag contentType
        );

        //- Add a (.vtp) file
        //  \return True if file is non-empty
        inline bool append_vtp(const fileName& file);

        //- Add a (.vtp) file with name
        //  \return True if file is non-empty
        inline bool append_vtp(const word& name, const fileName& file);

        //- Add a (.vtu) file
        //  \return True if file is non-empty
        inline bool append_vtu(const fileName& file);

        //- Add a (.vtu) file with name
        //  \return True if file is non-empty
        inline bool append_vtu(const word& name, const fileName& file);


    // Content Management

        //- Sanity fixes on the data
        void repair(bool collapse=false);

        //- Add in content from another vtm and place under the given block
        //- name.
        void add(const word& blockName, const vtmWriter& other);

        //- Add in content from another vtm and place under the given block
        //- name. Adjust the added 'file' entries to include the given prefix.
        void add
        (
            const word& blockName,
            const fileName& prefix,
            const vtmWriter& other
        );

    // Writing

        //- Open file for writing (creates parent directory) and write the
        //- blocks and TimeValue.
        //  The file name is with/without an extension.
        //  \return number of data sets
        label write(const fileName& file);

        //- Print debug view of block and dataset contents
        void dump(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtmWriterI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

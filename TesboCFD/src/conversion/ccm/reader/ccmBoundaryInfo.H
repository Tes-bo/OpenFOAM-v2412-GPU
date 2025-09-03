/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

Description
    Container for holding STARCCM boundary information

\*---------------------------------------------------------------------------*/

#ifndef ccmBoundaryInfo_H
#define ccmBoundaryInfo_H

#include "word.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//! \cond internal_use
// internal use only - skip doxygen documentation
namespace Foam
{
namespace ccm
{

class ccmBoundaryInfo;
Ostream& operator<<(Ostream& os, const ccmBoundaryInfo& entry);

/*---------------------------------------------------------------------------*\
                 Class Foam::ccm::ccmBoundaryInfo Declaration
\*---------------------------------------------------------------------------*/

//- Helper when reading raw boundary information
class ccmBoundaryInfo
{
public:
    //- The ccm region
    int ccmIndex;

    //- Number of faces
    label size;

    //- The openfoam patch id to map to
    label patchId;

    //- The patch name, as per the BoundaryRegion "Label" entry
    std::string patchName;

    //- The patch type, as per the BoundaryRegion "BoundaryType" entry
    std::string patchType;


    //- Construct null
    ccmBoundaryInfo()
    :
        ccmIndex(0),
        size(0),
        patchId(-1),
        patchName(),
        patchType("patch")
    {}


    //- Set patch name, default to "patch_CCMID" for an empty string
    void setPatchName(const std::string& str)
    {
        if (str.empty())
        {
            patchName = "patch_" + ::Foam::name(ccmIndex);
        }
        else
        {
            patchName = str;
        }
    }


    //- Info doesn't match if the ccm boundaries are different
    bool operator!=(const ccmBoundaryInfo& rhs) const
    {
        return ccmIndex != rhs.ccmIndex;
    }


    //- Ostream Operator
    friend Ostream& operator<<(Ostream& os, const ccmBoundaryInfo& entry)
    {
        os  << "BoundaryFaces-" << entry.ccmIndex;
        os  << " size=" << entry.size;
        os  << " name=" << entry.patchName;
        os  << " type=" << entry.patchType;
        os  << " foam-patch=" << entry.patchId;

        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ccm
} // End namespace Foam

//! \endcond
// internal use only - skip doxygen documentation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

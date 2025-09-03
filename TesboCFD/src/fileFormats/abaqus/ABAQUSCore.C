/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ABAQUSCore.H"
#include "IFstream.H"
#include "ListOps.H"
#include "stringOps.H"
#include "SpanStream.H"
#include "cellModel.H"
#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

static Foam::Map<Foam::labelList> abaqusToFoamFaceAddr_;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Use peek(), get(), unget() to detect and skip lines that appear to be
// "** comment" lines.
//
// Uses a mix of std::istream and ISstream methods
// since we need the low-level get()/unget()

static void skipComments(ISstream& iss)
{
    #if OPENFOAM < 2002
    string line;
    #endif

    auto& is = iss.stdStream();

    bool isComment = true;
    while (isComment)
    {
        isComment = ('*' == is.peek());

        if (isComment)
        {
            // Get and check the next one
            (void) is.get();

            isComment = ('*' == is.peek());
            if (isComment)
            {
                // Found "** ..." (a comment) - read/discard
                // ISstream::getLine to keep track of the line numbers

                #if OPENFOAM >= 2002
                iss.getLine(nullptr);
                #else
                iss.getLine(line);
                #endif
            }
            else
            {
                // Not a comment
                // - unget the '*', implicitly break out of loop
                is.unget();
            }
        }
    }

    iss.syncState();
}


// Get an identifier of the form "NSET=..." (case-insensitive)
// Return the string on success, an empty string on failure

static string getIdentifier(const word& keyword, string& inputLine)
{
    // Strip out whitespace (not a valid Abaqus identifier anyhow)
    // - makes parsing easier, avoids tab/carriage-returns etc.

    stringOps::inplaceRemoveSpace(inputLine);

    // Do string comparisons in upper-case

    const auto key(stringOps::upper(keyword));
    const auto line(stringOps::upper(inputLine));

    // Extract "..,key=value,key2=value,"

    // Not sure if we need the additional ',' prefix
    // in search to avoid similar keys.

    auto beg = line.find("," + key + "=");

    if (beg != std::string::npos)
    {
        // Skip past the '='
        beg += key.size() + 2;

        // The closing comma
        auto len = line.find(',', beg);
        if (len != std::string::npos)
        {
            len -= beg;
        }

        // Substring from inputLine (not uppercase!)
        return inputLine.substr(beg, len);
    }

    // Not found
    return string();
}


// Walk the string content (CSV format) to append integer labels
// until the line is exhausted or the list is full.
//
// Return false on read error or if the line exhausted while getting
// element.

static bool appendCsvLabels
(
    const std::string& line,
    labelUList& elemNodes,
    label& nodei
)
{
    const label nNodes = elemNodes.size();

    std::size_t pos = 0;

    while (nodei < nNodes && pos != std::string::npos)
    {
        auto beg = pos;
        auto len = line.find(',', pos);

        if (len == std::string::npos)
        {
            pos = len;
        }
        else
        {
            pos = len + 1;
            len -= beg;
        }

        if (readLabel(line.substr(beg, len), elemNodes[nodei]))
        {
            ++nodei;
        }
        else
        {
            // Read error, or need another line
            return false;
        }
    }

    return (nodei >= nNodes);
}


} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::Map<Foam::labelList>&
Foam::fileFormats::ABAQUSCore::abaqusToFoamFaceAddr()
{
    if (abaqusToFoamFaceAddr_.empty())
    {
        abaqusToFoamFaceAddr_.emplace(abaqusTet, labelList({3, 2, 0, 1}));
        abaqusToFoamFaceAddr_.emplace(abaqusPrism, labelList({0, 1, 4, 3, 2}));
        abaqusToFoamFaceAddr_.emplace(abaqusHex, labelList({4, 5, 2, 1, 3, 0}));
    }

    return abaqusToFoamFaceAddr_;
}


Foam::fileFormats::ABAQUSCore::shapeType
Foam::fileFormats::ABAQUSCore::getElementType(const string& elemTypeName)
{
    // Check for element-type.
    #undef  checkElemType
    #define checkElemType(Name) elemTypeName.contains(Name)

    if
    (
        checkElemType("S3")
     || checkElemType("CPE3")
     || checkElemType("2D3")
    )
    {
        return shapeType::abaqusTria;
    }
    else if
    (
        checkElemType("S4")
     || checkElemType("CPE4")
     || checkElemType("2D4")
     || checkElemType("CPEG4")
    )
    {
        return shapeType::abaqusQuad;
    }
    else if
    (
        checkElemType("3D4")    // C3D4*, Q3D4, ...
    )
    {
        return shapeType::abaqusTet;
    }
    else if
    (
        checkElemType("3D5")    // C3D5*
    )
    {
        return shapeType::abaqusPyr;
    }
    else if
    (
        checkElemType("3D6")    // C3D6*
    )
    {
        return shapeType::abaqusPrism;
    }
    else if
    (
        checkElemType("3D8")    // C3D8*
    )
    {
        return shapeType::abaqusHex;
    }

    #undef checkElemType

    return shapeType::abaqusUnknownShape;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label
Foam::fileFormats::ABAQUSCore::readHelper::addNewElset
(
    const std::string& setName
)
{
    if (elsetMap_.empty())
    {
        // Always have a lookup for empty string
        elsetMap_.set(string::null, 0);
    }

    if (setName.empty())
    {
        return 0;
    }

    // Direct case-sensitive lookup - it might be there

    label setId = elsetMap_.lookup(setName, -1);
    if (setId >= 0)
    {
        return setId;
    }


    // Case-insensitive search, use upper-case

    const auto needle(stringOps::upper(setName));

    forAllConstIters(elsetMap_, iter)
    {
        const auto haystack(stringOps::upper(iter.key()));

        if (needle == haystack)
        {
            return iter.val();
        }
    }

    // Not there. Save at the next location
    setId = elsetMap_.size();
    elsetMap_.set(setName, setId);

    return setId;
}


Foam::label
Foam::fileFormats::ABAQUSCore::readHelper::readPoints
(
    ISstream& is
)
{
    const label initialCount = points_.size();

    char sep; // Comma separator (dummy)
    string line;
    label id;
    point p;

    // Read nodes (points) until next "*Section"
    while (is.peek() != '*' && is.peek() != EOF)
    {
        // Grab the line and wrap as spanstream
        is.getLine(line);

        if (line.empty())
        {
            // Not sure if we should terminate on blank lines?
            continue;
        }
        ISpanStream ss(line);

        // Parse line for ID, X, Y, Z
        ss >> id >> sep >> p.x() >> sep >> p.y() >> sep >> p.z();

        nodeIds_.append(id);
        points_.append(p);
    }

    return (points_.size() - initialCount);
}


Foam::label
Foam::fileFormats::ABAQUSCore::readHelper::readElements
(
    ISstream& is,
    const ABAQUSCore::shapeType shape,
    const label setId
)
{
    // Info<< "*Element" << nl;

    const label nNodes = ABAQUSCore::nPoints(shape);

    if (!nNodes)
    {
        return 0;
    }

    const label initialCount = elemTypes_.size();

    char sep; // Comma separator (dummy)
    string line;
    label id;

    labelList elemNodes(nNodes, Zero);

    // Read element connectivity until next "*Section"

    // Parse for ID, node1, node2, ...
    while (is.peek() != '*' && is.peek() != EOF)
    {
        // elemNodes = Zero;   for sanity checks?

        is >> id >> sep;

        label nodei = 0;
        while (nodei < nNodes)
        {
            // Grab the rest of the line, walk through CSV fields
            is.getLine(line);

            appendCsvLabels(line, elemNodes, nodei);
        }

        // Checks?

        connectivity_.append(elemNodes);
        elemTypes_.append(shape);
        elemIds_.append(id);
        elsetIds_.append(setId);
    }

    return (elemTypes_.size() - initialCount);
}


Foam::label
Foam::fileFormats::ABAQUSCore::readHelper::readSurfaceElements
(
    ISstream& is,
    const label setId
)
{
    // Info<< "*Surface" << nl;

    // Models for supported solids (need to face mapping)
    const cellModel& tet   = cellModel::ref(cellModel::TET);
    const cellModel& prism = cellModel::ref(cellModel::PRISM);
    const cellModel& hex   = cellModel::ref(cellModel::HEX);

    // Face mapping from Abaqus cellModel to OpenFOAM cellModel
    const auto& abqToFoamFaceMap = abaqusToFoamFaceAddr();

    const label initialCount = elemTypes_.size();

    char sep; // Comma separator (dummy)
    string line;
    label id;

    // Read until next "*Section"

    // Parse for elemId, sideId.
    // Eg, "1235, S1"
    while (is.peek() != '*' && is.peek() != EOF)
    {
        is >> id >> sep;
        is.getLine(line);

        const word sideName(word::validate(stringOps::upper(line)));

        if
        (
            sideName.size() != 2
         || sideName[0] != 'S'
         || !std::isdigit(sideName[1])
        )
        {
            Info<< "Abaqus reader: unsupported surface element side "
                << id << ", " << sideName << nl;
            continue;
        }

        const label index = elemIds_.find(id);
        if (id <= 0 || index < 0)
        {
            Info<< "Abaqus reader: unsupported surface element "
                << id << nl;
            continue;
        }

        const auto faceIdIter = abqToFoamFaceMap.cfind(elemTypes_[index]);
        if (!faceIdIter.good())
        {
            Info<< "Abaqus reader: reject non-solid shape: " << nl;
        }

        // The abaqus element side number (1-based)
        const label sideNum = (sideName[1] - '0');

        const label foamFaceNum = (*faceIdIter)[sideNum - 1];

        const labelList& connect = connectivity_[index];

        // Nodes for the derived shell element
        labelList elemNodes;

        switch (elemTypes_[index])
        {
            case shapeType::abaqusTet:
            {
                elemNodes = labelList(connect, tet.modelFaces()[foamFaceNum]);
                break;
            }
            case shapeType::abaqusPrism:
            {
                elemNodes = labelList(connect, prism.modelFaces()[foamFaceNum]);
                break;
            }
            case shapeType::abaqusHex:
            {
                elemNodes = labelList(connect, hex.modelFaces()[foamFaceNum]);
                break;
            }
            default:
                break;
        }

        enum shapeType shape = shapeType::abaqusUnknownShape;

        if (elemNodes.size() == 3)
        {
            shape = shapeType::abaqusTria;
        }
        else if (elemNodes.size() == 4)
        {
            shape = shapeType::abaqusQuad;
        }
        else
        {
            // Cannot happen
            FatalErrorInFunction
                << "Could not map face side for "
                << id << ", " << sideName << nl
                << exit(FatalError);
        }

        // Synthesize face Id from solid element Id and side Id
        const label newElemId = ABAQUSCore::encodeSolidId(id, sideNum);

        // Further checks?
        connectivity_.append(std::move(elemNodes));
        elemTypes_.append(shape);
        elemIds_.append(newElemId);
        elsetIds_.append(setId);
    }

    return (elemTypes_.size() - initialCount);
}


void Foam::fileFormats::ABAQUSCore::readHelper::read
(
    ISstream& is
)
{
    clear();

    label nread;
    string line;

    while (is.good())
    {
        is.getLine(line);

        // Start processing on "*Section-Name",
        // but skip "** comments" etc
        if (line[0] != '*' || !std::isalpha(line[1]))
        {
            continue;
        }

        // Some abaqus files use upper-case or mixed-case for section names,
        // convert all to upper-case for ease.

        const string upperLine(stringOps::upper(line));

        //
        // "*Nodes" section
        //
        if (upperLine.starts_with("*NODE"))
        {
            // Ignore "NSET=...", we cannot do anything useful with it

            skipComments(is);

            nread = readPoints(is);

            if (verbose_)
            {
                InfoErr
                    << "Read " << nread << " *NODE entries" << nl;
            }
            continue;
        }

        //
        // "*Element" section
        //
        if (upperLine.starts_with("*ELEMENT,"))
        {
            // Must have "TYPE=..."
            const string elemTypeName(getIdentifier("TYPE", line));

            // May have "ELSET=..." on the same line
            const string elsetName(getIdentifier("ELSET", line));

            const shapeType shape(getElementType(elemTypeName));

            if (!ABAQUSCore::nPoints(shape))
            {
                // Unknown/unsupported
                if (verbose_)
                {
                    InfoErr
                        << "Ignore abaqus element type: "
                        << elemTypeName << nl;
                }
                continue;
            }

            const label elsetId = addNewElset(elsetName);

            skipComments(is);

            nread = readElements(is, shape, elsetId);

            if (verbose_)
            {
                InfoErr
                    << "Read " << nread << " *ELEMENT entries ("
                    << elemTypeName << ") elset="
                    << elsetName << nl;
            }
            continue;
        }

        //
        // "*Surface" section
        //
        if (upperLine.starts_with("*SURFACE,"))
        {
            // Require "NAME=..." on the same line
            const string elsetName(getIdentifier("NAME", line));

            // May have "TYPE=..." on the same line.
            // If missing, default is ELEMENT.
            const string surfTypeName(getIdentifier("TYPE", line));

            if
            (
                !surfTypeName.empty()
             && stringOps::upper(surfTypeName) != "ELEMENT"
            )
            {
                Info<< "Reading abaqus surface type "
                    << surfTypeName << " is not implemented" << nl;
                continue;
            }

            // Treat like an element set
            const label elsetId = addNewElset(elsetName);

            skipComments(is);

            nread = readSurfaceElements(is, elsetId);

            if (verbose_)
            {
                InfoErr
                    << "Read " << nread << " *SURFACE entries for "
                    << elsetName << nl;
            }
            continue;
        }
    }
}


void Foam::fileFormats::ABAQUSCore::readHelper::purge_solids()
{
    // Negative set
    bitSet select(elemTypes_.size(), false);

    forAll(elemTypes_, i)
    {
        if (!isValidType(elemTypes_[i]) || isSolidType(elemTypes_[i]))
        {
            select.set(i);
        }
    }

    if (select.any())
    {
        select.flip();

        inplaceSubset(select, connectivity_);
        inplaceSubset(select, elemTypes_);

        inplaceSubset(select, elemIds_);
        inplaceSubset(select, elsetIds_);
    }
}


void Foam::fileFormats::ABAQUSCore::readHelper::compact_nodes()
{
    if (!nodeIds_.empty())
    {
        // Has original 1-based ids
        //
        // Need to convert to local (0-based) points
        // in the order in which we read them
        // and compact unused values

        // Could construct a sort order to preserve the original
        // point order, but that is not likely relevant for anyone.


        // Which original node ids actually being used by elements?

        // We may have many ids, but speculate that they are sparse
        // and have high element numbers.

        // Use a Map instead of labelList.

        Map<label> nodeIdRemapping(2*points_.size());

        // Pass 1: which nodes are being used?

        for (const labelList& elem : connectivity_)
        {
            for (const label origId : elem)
            {
                nodeIdRemapping(origId) = 0; // any value
            }
        }

        // Define compact local points, finalize the node id remapping

        label nPoints = 0;
        labelList oldToNewLocal(nodeIds_.size(), -1);

        forAll(nodeIds_, i)
        {
            const label origId = nodeIds_[i];

            if (nodeIdRemapping.found(origId))
            {
                oldToNewLocal[i] = nPoints;
                nodeIdRemapping(origId) = nPoints;
                ++nPoints;
            }
        }

        // Prune out -1 values (shrinks list)
        inplaceReorder(oldToNewLocal, points_, true);

        // Relabel the elements

        for (labelList& elem : connectivity_)
        {
            for (label& id : elem)
            {
                id = nodeIdRemapping[id];
            }
        }

        // Done!
        nodeIds_.clear();
    }
    else
    {
        // Already numbered (0-based), but perhaps not compacted

        // Which node ids actually being used by elements?
        bitSet usedNodeIds(points_.size());

        for (const labelList& elem : connectivity_)
        {
            usedNodeIds.set(elem);
        }

        // Compact the numbers
        labelList oldToNewLocal = invert(points_.size(), usedNodeIds);

        // Prune out -1 values (shrinks list)
        inplaceReorder(oldToNewLocal, points_, true);


        // Renumber non-compact to compact
        for (labelList& elem : connectivity_)
        {
            inplaceRenumber(oldToNewLocal, elem);
        }
    }
}


void Foam::fileFormats::ABAQUSCore::readHelper::renumber_elements_1to0()
{
    for (label& elemId : elemIds_)
    {
        renumber0_elemId(elemId);
    }
}


void Foam::fileFormats::ABAQUSCore::writePoints
(
    Ostream& os,
    const UList<point>& points,
    const scalar scaleFactor
)
{
    if (points.empty())
    {
        return;
    }

    // Set the precision of the points data to 10
    os.precision(10);

    // Force decimal point for Fortran input
    os.setf(std::ios::showpoint);

    label vertId = 1;  // 1-based vertex labels

    os << "*NODE" << nl;
    for (const point& p : points)
    {
        // Convert [m] -> [mm] etc
        os  << "  "
            << vertId << ", "
            << (scaleFactor * p.x()) << ','
            << (scaleFactor * p.y()) << ','
            << (scaleFactor * p.z()) << nl;

        ++vertId;
    }
}


Foam::label Foam::fileFormats::ABAQUSCore::faceDecomposition
(
    const UList<point>& points,
    const UList<face>& faces,
    labelList& decompOffsets,
    DynamicList<face>& decompFaces
)
{
    // On-demand face decomposition (triangulation)

    decompOffsets.resize(faces.size()+1);
    decompFaces.clear();

    auto offsetIter = decompOffsets.begin();
    *offsetIter = 0; // The first offset is always zero

    for (const face& f : faces)
    {
        const label n = f.size();

        if (n != 3 && n != 4)
        {
            // Decompose non-tri/quad into tris
            f.triangles(points, decompFaces);
        }

        // The end offset, which is the next begin offset
        *(++offsetIter) = decompFaces.size();
    }

    return decompFaces.size();
}


// ************************************************************************* //

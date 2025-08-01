// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkEnSightGoldBinaryReader.h"

#include "vtkByteSwap.h"
#include "vtkCellData.h"
#include "vtkCharArray.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkImageData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtksys/Encoding.hxx"
#include "vtksys/FStream.hxx"
#include "vtksys/RegularExpression.hxx"
#include "vtksys/SystemTools.hxx"

#include <algorithm> /* std::remove */
#include <array>
#include <cctype>
#include <map>
#include <numeric>
#include <string>
#include <sys/stat.h>
#include <vector>

#if defined(_WIN32)
#define VTK_STAT_STRUCT struct _stat64
#define VTK_STAT_FUNC _stat64
#elif defined _DARWIN_FEATURE_64_BIT_INODE || defined __FreeBSD__ || defined __NetBSD__ ||         \
  defined __OpenBSD__
// The BSDs use stat().
#define VTK_STAT_STRUCT struct stat
#define VTK_STAT_FUNC stat
#elif defined __EMSCRIPTEN__
#if defined _LARGEFILE64_SOURCE
#define VTK_STAT_STRUCT struct stat64
#define VTK_STAT_FUNC stat64
#else
#define VTK_STAT_STRUCT struct stat
#define VTK_STAT_FUNC stat
#endif
#else
// here, we're relying on _FILE_OFFSET_BITS defined in vtkWin32Header.h to help
// us on POSIX without resorting to using stat64.
#define VTK_STAT_STRUCT struct stat64
#define VTK_STAT_FUNC stat64
#endif

VTK_ABI_NAMESPACE_BEGIN
class vtkEnSightGoldBinaryReader::vtkUtilities
{
  static int GetDestinationComponent(int srcComponent, int numComponents)
  {
    if (numComponents == 6)
    {
      // for 6 component tensors, the symmetric tensor components XZ and YZ are interchanged
      // see #10637.
      switch (srcComponent)
      {
        case 4:
          return 5;

        case 5:
          return 4;
      }
    }

    return srcComponent;
  }

public:
  static vtkSmartPointer<vtkFloatArray> ReadVariableFloats(const char* sectionHeader,
    vtkEnSightGoldBinaryReader* self, const char* description, vtkDataSetAttributes* dsa,
    vtkIdType numElements, int numComponents, int component = -1)
  {
    vtksys::RegularExpression regEx("^[^ ]+ ([^ ]+)");
    const bool match = regEx.find(sectionHeader);
    const bool hasUndef = (match && regEx.match(1) == "undef");
    const bool hasPartial = (match && regEx.match(1) == "partial");

    float undefValue{ 0 };
    if (hasUndef)
    {
      self->ReadFloat(&undefValue);
    }

    vtkNew<vtkIdList> partialIndices;
    if (hasPartial)
    {
      int count;
      self->ReadInt(&count);

      std::vector<int> buffer(count);
      self->ReadIntArray(&buffer.front(), count);

      partialIndices->SetNumberOfIds(count);
      std::transform(
        buffer.begin(), buffer.end(), partialIndices->GetPointer(0), [](vtkIdType val) {
          return val - 1; /* since ensight indices start with 1*/
        });
    }

    // replace undefined values with "internal undef" which in ParaView is NaN
    auto replaceUndef = [&](vtkFloatArray* farray) {
      if (hasUndef)
      {
        const float nanfloat = std::nanf("1");
        for (vtkIdType cc = 0; cc < numElements; ++cc)
        {
          if (farray->GetTypedComponent(cc, 0) == undefValue)
          {
            farray->SetTypedComponent(cc, 0, nanfloat);
          }
        }
      }
    };

    auto readComponent = [&](vtkIdType count) {
      vtkNew<vtkFloatArray> buffer;
      buffer->SetNumberOfTuples(count);
      if (hasPartial)
      {
        // fill with NaNs
        buffer->FillValue(std::nanf("1"));

        vtkNew<vtkFloatArray> pbuffer;
        pbuffer->SetNumberOfTuples(partialIndices->GetNumberOfIds());
        self->ReadFloatArray(
          pbuffer->GetPointer(0), static_cast<int>(partialIndices->GetNumberOfIds()));

        // now copy the tuples over from pbuffer to buffer.
        vtkNew<vtkIdList> srcIds;
        srcIds->SetNumberOfIds(partialIndices->GetNumberOfIds());
        std::iota(srcIds->begin(), srcIds->end(), 0);
        buffer->InsertTuples(partialIndices, srcIds, pbuffer);
      }
      else
      {
        self->ReadFloatArray(buffer->GetPointer(0), count);
        replaceUndef(buffer);
      }
      return buffer;
    };

    if (numComponents == 1)
    {
      return readComponent(numElements);
    }
    else if (numComponents > 1 && component != -1)
    {
      vtkSmartPointer<vtkFloatArray> array;
      if (component == 0)
      {
        array = vtk::TakeSmartPointer(vtkFloatArray::New());
        array->SetNumberOfComponents(numComponents);
        array->SetNumberOfTuples(numElements);
      }
      else
      {
        array = vtkFloatArray::SafeDownCast(dsa->GetArray(description));
        assert(array && array->GetNumberOfComponents() == numComponents);
      }
      auto buffer = readComponent(numElements);
      array->CopyComponent(component, buffer, 0);
      return array;
    }
    else if (numComponents > 1 && component == -1)
    {
      vtkSmartPointer<vtkFloatArray> array;
      array = vtk::TakeSmartPointer(vtkFloatArray::New());
      array->SetNumberOfComponents(numComponents);
      array->SetNumberOfTuples(numElements);
      for (int comp = 0; comp < numComponents; ++comp)
      {
        const int destComponent = vtkUtilities::GetDestinationComponent(comp, numComponents);
        auto buffer = readComponent(numElements);
        array->CopyComponent(destComponent, buffer, 0);
      }
      return array;
    }

    return nullptr;
  }
};

vtkStandardNewMacro(vtkEnSightGoldBinaryReader);
class vtkEnSightGoldBinaryReader::FileOffsetMapInternal
{
  typedef std::string MapKey;
  typedef std::map<int, vtkTypeInt64> MapValue;

public:
  typedef std::map<MapKey, MapValue>::const_iterator const_iterator;
  typedef std::map<MapKey, MapValue>::value_type value_type;

  std::map<MapKey, MapValue> Map;
};

// This is half the precision of an int.
#define MAXIMUM_PART_ID 65536

//------------------------------------------------------------------------------
vtkEnSightGoldBinaryReader::vtkEnSightGoldBinaryReader()
{
  this->FileOffsets = new vtkEnSightGoldBinaryReader::FileOffsetMapInternal;

  this->GoldIFile = nullptr;
  this->FileSize = 0;
  this->SizeOfInt = sizeof(int);
  this->Fortran = 0;
  this->FortranSkipBytes = 0;
  this->NodeIdsListed = 0;
  this->ElementIdsListed = 0;
}

//------------------------------------------------------------------------------
vtkEnSightGoldBinaryReader::~vtkEnSightGoldBinaryReader()
{
  delete this->FileOffsets;
  delete this->GoldIFile;
  this->GoldIFile = nullptr;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::OpenFile(const char* filename)
{
  if (!filename)
  {
    vtkErrorMacro(<< "Missing filename.");
    return 0;
  }

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  // Open the new file
  vtkDebugMacro(<< "Opening file " << filename);
  VTK_STAT_STRUCT fs;
  if (!VTK_STAT_FUNC(filename, &fs))
  {
    // Find out how big the file is.
    this->FileSize = static_cast<vtkTypeUInt64>(fs.st_size);

    std::ios_base::openmode mode = ios::in;
#ifdef _WIN32
    mode |= ios::binary;
#endif
    this->GoldIFile = new vtksys::ifstream(filename, mode);
  }
  else
  {
    vtkErrorMacro("stat failed.");
    return 0;
  }
  if (!this->GoldIFile || this->GoldIFile->fail())
  {
    vtkErrorMacro(<< "Could not open file " << filename);
    return 0;
  }

  // we now need to check for Fortran and byte ordering

  // we need to look at the first 4 bytes of the file, and the 84-87 bytes
  // of the file to correctly determine what it is. If we only check the first
  // 4 bytes we can get incorrect detection if it is a property file named "P"
  // we check the 84-87 bytes as that is the start of the next line on a fortran file

  char result[88];
  this->GoldIFile->read(result, 88);
  if (this->GoldIFile->eof() || this->GoldIFile->fail())
  {
    vtkErrorMacro(<< filename << " is missing header information");
    return 0;
  }
  this->GoldIFile->seekg(0, ios::beg); // reset the file to the start

  // if the first 4 bytes is the length, then this data is no doubt
  // a fortran data write!, copy the last 76 into the beginning
  char le_len[4] = { 0x50, 0x00, 0x00, 0x00 };
  char be_len[4] = { 0x00, 0x00, 0x00, 0x50 };

  // the fortran test here depends on the byte ordering. But if the user didn't
  // set any byte ordering then, we have to try both byte orderings. There was a
  // bug here which was resulting in binary-fortran-big-endian files being read
  // incorrectly on intel machines (BUG #10593). This dual-check avoids that
  // bug.
  bool le_isFortran = true;
  bool be_isFortran = true;
  for (int c = 0; c < 4; c++)
  {
    le_isFortran = le_isFortran && (result[c] == le_len[c]) && (result[c + 84] == le_len[c]);
    be_isFortran = be_isFortran && (result[c] == be_len[c]) && (result[c + 84] == be_len[c]);
  }

  switch (this->ByteOrder)
  {
    case FILE_BIG_ENDIAN:
      this->Fortran = be_isFortran;
      break;

    case FILE_LITTLE_ENDIAN:
      this->Fortran = le_isFortran;
      break;

    case FILE_UNKNOWN_ENDIAN:
      if (le_isFortran)
      {
        this->Fortran = true;
        this->ByteOrder = FILE_LITTLE_ENDIAN;
      }
      else if (be_isFortran)
      {
        this->Fortran = true;
        this->ByteOrder = FILE_BIG_ENDIAN;
      }
      else
      {
        this->Fortran = false;
      }
      break;
  }
  this->FortranSkipBytes = this->Fortran ? 8 : 0;

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::InitializeFile(const char* fileName)
{
  char line[80], subLine[80];

  // Initialize
  //
  if (!fileName)
  {
    vtkErrorMacro("A GeometryFileName must be specified in the case file.");
    return 0;
  }
  std::string sfilename;
  std::string filenameString(fileName);
  this->SanitizeFileName(filenameString);
  if (this->FilePath)
  {
    sfilename = this->FilePath;
    if (sfilename.at(sfilename.length() - 1) != '/')
    {
      sfilename += "/";
    }
    sfilename += filenameString;
    vtkDebugMacro("full path to geometry file: " << sfilename);
  }
  else
  {
    sfilename = filenameString;
  }

  if (this->OpenFile(sfilename.c_str()) == 0)
  {
    vtkErrorMacro("Unable to open file: " << sfilename);
    return 0;
  }

  line[0] = '\0';
  subLine[0] = '\0';
  if (this->ReadLine(line) == 0)
  {
    vtkErrorMacro("Error with line reading upon file initialization");
    return 0;
  }

  if (sscanf(line, " %*s %s", subLine) != 1)
  {
    vtkErrorMacro("Error with subline extraction upon file initialization");
    return 0;
  }

  if (strncmp(subLine, "Binary", 6) != 0 && strncmp(subLine, "binary", 6) != 0)
  {
    vtkErrorMacro("This is not a binary data set. Try "
      << "vtkEnSightGoldReader.");
    return 0;
  }
  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadGeometryFile(
  const char* fileName, int timeStep, vtkMultiBlockDataSet* output)
{
  char line[80], subLine[80], nameline[80];
  int partId, realId;
  int lineRead, i;

  if (!this->InitializeFile(fileName))
  {
    return 0;
  }

  // this will close the file, so we need to reinitialize it
  int numberOfTimeStepsInFile = this->CountTimeSteps();

  if (!this->InitializeFile(fileName))
  {
    return 0;
  }

  if (this->UseFileSets)
  {
    if (numberOfTimeStepsInFile > 1)
    {
      this->AddFileIndexToCache(fileName);

      i = this->SeekToCachedTimeStep(fileName, timeStep - 1);
      // start w/ the number of TS we skipped, not the one we are at
      // if we are not at the appropriate time step yet, we keep searching
      for (; i < timeStep - 1; i++)
      {
        if (!this->SkipTimeStep())
        {
          return 0;
        }
      }
    }

    // use do-while here to initialize 'line' before 'strncmp' is applied
    // Thanks go to Brancois for care of this issue
    do
    {
      this->ReadLine(line);
    } while (strncmp(line, "BEGIN TIME STEP", 15) != 0);
    // found a time step -> cache it
    this->AddTimeStepToCache(fileName, timeStep - 1, this->GoldIFile->tellg());
  }

  // Skip the 2 description lines.
  this->ReadLine(line);
  this->ReadLine(line);

  // Read the node id and element id lines.
  this->ReadLine(line);
  sscanf(line, " %*s %*s %s", subLine);
  if (strncmp(subLine, "given", 5) == 0)
  {
    this->NodeIdsListed = 1;
  }
  else if (strncmp(subLine, "ignore", 6) == 0)
  {
    this->NodeIdsListed = 1;
  }
  else
  {
    this->NodeIdsListed = 0;
  }

  this->ReadLine(line);
  sscanf(line, " %*s %*s %s", subLine);
  if (strncmp(subLine, "given", 5) == 0)
  {
    this->ElementIdsListed = 1;
  }
  else if (strncmp(subLine, "ignore", 6) == 0)
  {
    this->ElementIdsListed = 1;
  }
  else
  {
    this->ElementIdsListed = 0;
  }

  lineRead = this->ReadLine(line); // "extents" or "part"
  if (strncmp(line, "extents", 7) == 0)
  {
    // Skipping the extents.
    this->GoldIFile->seekg(6 * sizeof(float) + this->FortranSkipBytes, ios::cur);
    lineRead = this->ReadLine(line); // "part"
  }

  while (lineRead > 0 && strncmp(line, "part", 4) == 0)
  {
    this->ReadPartId(&partId);
    partId--; // EnSight starts #ing at 1.
    if (partId < 0 || partId >= MAXIMUM_PART_ID)
    {
      vtkErrorMacro("Invalid part id; check that ByteOrder is set correctly.");
      return 0;
    }
    realId = this->InsertNewPartId(partId);

    // Increment the number of geometry parts such that the measured geometry,
    // if any, can be properly combined into a vtkMultiBlockDataSet object.
    // --- fix to bug #7453
    this->NumberOfGeometryParts++;

    this->ReadLine(line); // part description line

    strncpy(nameline, line, 80); // 80 characters in line are allowed
    nameline[79] = '\0';         // Ensure nullptr character at end of part name
    char* name = strdup(nameline);

    // fix to bug #0008237
    // The original "return 1" operation upon "strncmp(line, "interface", 9) == 0"
    // was removed here as 'interface' is NOT a keyword of an EnSight Gold file.

    this->ReadLine(line);

    if (strncmp(line, "block", 5) == 0)
    {
      if (sscanf(line, " %*s %s", subLine) == 1)
      {
        if (strncmp(subLine, "rectilinear", 11) == 0)
        {
          // block rectilinear
          lineRead = this->CreateRectilinearGridOutput(realId, line, name, output);
        }
        else if (strncmp(subLine, "uniform", 7) == 0)
        {
          // block uniform
          lineRead = this->CreateImageDataOutput(realId, line, name, output);
        }
        else
        {
          // block iblanked
          lineRead = this->CreateStructuredGridOutput(realId, line, name, output);
        }
      }
      else
      {
        // block
        lineRead = this->CreateStructuredGridOutput(realId, line, name, output);
      }
    }
    else
    {
      lineRead = this->CreateUnstructuredGridOutput(realId, line, name, output);
      if (lineRead < 0)
      {
        free(name);
        delete this->GoldIFile;
        this->GoldIFile = nullptr;
        return 0;
      }
    }
    free(name);
  }

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  if (lineRead < 0)
  {
    return 0;
  }

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::CountTimeSteps()
{
  int count = 0;
  while (true)
  {
    int result = this->SkipTimeStep();
    if (result)
    {
      count++;
    }
    else
    {
      break;
    }
  }
  return count;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::SkipTimeStep()
{
  char line[80], subLine[80];
  int lineRead;

  line[0] = '\0';
  while (strncmp(line, "BEGIN TIME STEP", 15) != 0)
  {
    if (!this->ReadLine(line))
    {
      return 0;
    }
  }

  // Skip the 2 description lines.
  this->ReadLine(line);
  this->ReadLine(line);

  // Read the node id and element id lines.
  this->ReadLine(line);
  sscanf(line, " %*s %*s %s", subLine);
  if (strncmp(subLine, "given", 5) == 0 || strncmp(subLine, "ignore", 6) == 0)
  {
    this->NodeIdsListed = 1;
  }
  else
  {
    this->NodeIdsListed = 0;
  }

  this->ReadLine(line);
  sscanf(line, " %*s %*s %s", subLine);
  if (strncmp(subLine, "given", 5) == 0)
  {
    this->ElementIdsListed = 1;
  }
  else if (strncmp(subLine, "ignore", 6) == 0)
  {
    this->ElementIdsListed = 1;
  }
  else
  {
    this->ElementIdsListed = 0;
  }

  lineRead = this->ReadLine(line); // "extents" or "part"
  if (strncmp(line, "extents", 7) == 0)
  {
    // Skipping the extents.
    this->GoldIFile->seekg(6 * sizeof(float) + this->FortranSkipBytes, ios::cur);
    lineRead = this->ReadLine(line); // "part"
  }

  while (lineRead > 0 && strncmp(line, "part", 4) == 0)
  {
    int tmpInt;
    this->ReadPartId(&tmpInt);
    if (tmpInt < 0 || tmpInt > MAXIMUM_PART_ID)
    {
      vtkErrorMacro("Invalid part id; check that ByteOrder is set correctly.");
      return 0;
    }
    this->ReadLine(line); // part description line
    this->ReadLine(line);

    if (strncmp(line, "block", 5) == 0)
    {
      if (sscanf(line, " %*s %s", subLine) == 1)
      {
        if (strncmp(subLine, "rectilinear", 11) == 0)
        {
          // block rectilinear
          lineRead = this->SkipRectilinearGrid(line);
        }
        else if (strncmp(subLine, "uniform,", 7) == 0)
        {
          // block uniform
          lineRead = this->SkipImageData(line);
        }
        else
        {
          // block iblanked
          lineRead = this->SkipStructuredGrid(line);
        }
      }
      else
      {
        // block
        lineRead = this->SkipStructuredGrid(line);
      }
    }
    else
    {
      lineRead = this->SkipUnstructuredGrid(line);
    }
  }

  if (lineRead < 0)
  {
    delete this->GoldIFile;
    this->GoldIFile = nullptr;
    return 0;
  }

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::SkipStructuredGrid(char line[256])
{
  char subLine[80];
  int lineRead;
  int iblanked = 0;
  int dimensions[3];
  int numPts;

  if (sscanf(line, " %*s %s", subLine) == 1)
  {
    if (strncmp(subLine, "iblanked", 8) == 0)
    {
      iblanked = 1;
    }
  }

  this->ReadIntArray(dimensions, 3);
  numPts = dimensions[0] * dimensions[1] * dimensions[2];
  if (dimensions[0] < 0 ||
    static_cast<unsigned int>(dimensions[0] * this->SizeOfInt) > this->FileSize ||
    dimensions[1] < 0 ||
    static_cast<unsigned int>(dimensions[1] * this->SizeOfInt) > this->FileSize ||
    dimensions[2] < 0 ||
    static_cast<unsigned int>(dimensions[2] * this->SizeOfInt) > this->FileSize || numPts < 0 ||
    static_cast<unsigned int>(numPts * this->SizeOfInt) > this->FileSize)
  {
    vtkErrorMacro("Invalid dimensions read; check that ByteOrder is set correctly.");
    return -1;
  }

  // Skip xCoords, yCoords and zCoords.
  this->GoldIFile->seekg(sizeof(float) * numPts * 3 + this->FortranSkipBytes, ios::cur);

  if (iblanked)
  { // skip iblank array.
    this->GoldIFile->seekg(numPts * sizeof(int) + this->FortranSkipBytes, ios::cur);
  }

  // reading next line to check for EOF
  lineRead = this->ReadLine(line);
  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::SkipUnstructuredGrid(char line[256])
{
  int lineRead = 1;
  int i;
  int numElements;
  int cellType;

  while (lineRead && strncmp(line, "part", 4) != 0)
  {
    if (strncmp(line, "coordinates", 11) == 0)
    {
      vtkDebugMacro("coordinates");
      int numPts;

      this->ReadInt(&numPts);
      if (numPts < 0 || static_cast<unsigned int>(numPts * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of points; check that ByteOrder is set correctly.");
        return -1;
      }

      vtkDebugMacro("num. points: " << numPts);

      if (this->NodeIdsListed)
      { // skip node ids.
        this->GoldIFile->seekg(sizeof(int) * numPts + this->FortranSkipBytes, ios::cur);
      }

      // Skip xCoords, yCoords and zCoords.
      this->GoldIFile->seekg(sizeof(float) * 3 * numPts + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "point", 5) == 0 || strncmp(line, "g_point", 7) == 0)
    {
      vtkDebugMacro("point");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of point cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      { // skip element ids.
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "bar2", 4) == 0 || strncmp(line, "g_bar2", 6) == 0)
    {
      vtkDebugMacro("bar2");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of bar2 cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * 2 * numElements + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "bar3", 4) == 0 || strncmp(line, "g_bar3", 6) == 0)
    {
      vtkDebugMacro("bar3");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of bar3 cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * 3 * numElements + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "nsided", 6) == 0 || strncmp(line, "g_nsided", 8) == 0)
    {
      vtkDebugMacro("nsided");
      int* numNodesPerElement;
      int numNodes = 0;

      // cellType = vtkEnSightReader::NSIDED;
      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of nsided cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      numNodesPerElement = new int[numElements];
      this->ReadIntArray(numNodesPerElement, numElements);
      for (i = 0; i < numElements; i++)
      {
        numNodes += numNodesPerElement[i];
      }
      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * numNodes + this->FortranSkipBytes, ios::cur);
      delete[] numNodesPerElement;
    }
    else if (strncmp(line, "tria3", 5) == 0 || strncmp(line, "tria6", 5) == 0 ||
      strncmp(line, "g_tria3", 7) == 0 || strncmp(line, "g_tria6", 7) == 0)
    {
      if (strncmp(line, "tria6", 5) == 0 || strncmp(line, "g_tria6", 7) == 0)
      {
        vtkDebugMacro("tria6");
        cellType = vtkEnSightReader::TRIA6;
      }
      else
      {
        vtkDebugMacro("tria3");
        cellType = vtkEnSightReader::TRIA3;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of triangle cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::TRIA6)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 6 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 3 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "quad4", 5) == 0 || strncmp(line, "quad8", 5) == 0 ||
      strncmp(line, "g_quad4", 7) == 0 || strncmp(line, "g_quad8", 7) == 0)
    {
      if (strncmp(line, "quad8", 5) == 0 || strncmp(line, "g_quad8", 7) == 0)
      {
        vtkDebugMacro("quad8");
        cellType = vtkEnSightReader::QUAD8;
      }
      else
      {
        vtkDebugMacro("quad4");
        cellType = vtkEnSightReader::QUAD4;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of quad cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::QUAD8)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 8 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 4 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "nfaced", 6) == 0)
    {
      vtkDebugMacro("nfaced");
      int* numFacesPerElement;
      int* numNodesPerFace;
      int numFaces = 0;
      int numNodes = 0;

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of nfaced cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      numFacesPerElement = new int[numElements];
      this->ReadIntArray(numFacesPerElement, numElements);
      for (i = 0; i < numElements; i++)
      {
        numFaces += numFacesPerElement[i];
      }
      delete[] numFacesPerElement;
      numNodesPerFace = new int[numFaces];
      this->ReadIntArray(numNodesPerFace, numFaces);
      for (i = 0; i < numFaces; i++)
      {
        numNodes += numNodesPerFace[i];
      }
      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * numNodes + this->FortranSkipBytes, ios::cur);
      delete[] numNodesPerFace;
    }
    else if (strncmp(line, "tetra4", 6) == 0 || strncmp(line, "tetra10", 7) == 0 ||
      strncmp(line, "g_tetra4", 8) == 0 || strncmp(line, "g_tetra10", 9) == 0)
    {
      if (strncmp(line, "tetra10", 7) == 0 || strncmp(line, "g_tetra10", 9) == 0)
      {
        vtkDebugMacro("tetra10");
        cellType = vtkEnSightReader::TETRA10;
      }
      else
      {
        vtkDebugMacro("tetra4");
        cellType = vtkEnSightReader::TETRA4;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro(
          "Invalid number of tetrahedral cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::TETRA10)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 10 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 4 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "pyramid5", 8) == 0 || strncmp(line, "pyramid13", 9) == 0 ||
      strncmp(line, "g_pyramid5", 10) == 0 || strncmp(line, "g_pyramid13", 11) == 0)
    {
      if (strncmp(line, "pyramid13", 9) == 0 || strncmp(line, "g_pyramid13", 11) == 0)
      {
        vtkDebugMacro("pyramid13");
        cellType = vtkEnSightReader::PYRAMID13;
      }
      else
      {
        vtkDebugMacro("pyramid5");
        cellType = vtkEnSightReader::PYRAMID5;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of pyramid cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::PYRAMID13)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 13 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 5 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "hexa8", 5) == 0 || strncmp(line, "hexa20", 6) == 0 ||
      strncmp(line, "g_hexa8", 7) == 0 || strncmp(line, "g_hexa20", 8) == 0)
    {
      if (strncmp(line, "hexa20", 6) == 0 || strncmp(line, "g_hexa20", 8) == 0)
      {
        vtkDebugMacro("hexa20");
        cellType = vtkEnSightReader::HEXA20;
      }
      else
      {
        vtkDebugMacro("hexa8");
        cellType = vtkEnSightReader::HEXA8;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of hexahedral cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::HEXA20)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 20 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 8 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "penta6", 6) == 0 || strncmp(line, "penta15", 7) == 0 ||
      strncmp(line, "g_penta6", 8) == 0 || strncmp(line, "g_penta15", 9) == 0)
    {
      if (strncmp(line, "penta15", 7) == 0 || strncmp(line, "g_penta15", 9) == 0)
      {
        vtkDebugMacro("penta15");
        cellType = vtkEnSightReader::PENTA15;
      }
      else
      {
        vtkDebugMacro("penta6");
        cellType = vtkEnSightReader::PENTA6;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of pentagonal cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::PENTA15)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 15 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 6 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "END TIME STEP", 13) == 0)
    {
      return 1;
    }
    else
    {
      vtkErrorMacro("undefined geometry file line");
      return -1;
    }
    lineRead = this->ReadLine(line);
  }
  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::SkipRectilinearGrid(char line[256])
{
  char subLine[80];
  int lineRead;
  int iblanked = 0;
  int dimensions[3];
  int numPts;

  if (sscanf(line, " %*s %*s %s", subLine) == 1)
  {
    if (strncmp(subLine, "iblanked", 8) == 0)
    {
      iblanked = 1;
    }
  }

  this->ReadIntArray(dimensions, 3);
  if (dimensions[0] < 0 ||
    static_cast<unsigned int>(dimensions[0] * this->SizeOfInt) > this->FileSize ||
    dimensions[1] < 0 ||
    static_cast<unsigned int>(dimensions[1] * this->SizeOfInt) > this->FileSize ||
    dimensions[2] < 0 ||
    static_cast<unsigned int>(dimensions[2] * this->SizeOfInt) > this->FileSize ||
    (dimensions[0] + dimensions[1] + dimensions[2]) < 0 ||
    static_cast<unsigned int>((dimensions[0] + dimensions[1] + dimensions[2]) * this->SizeOfInt) >
      this->FileSize)
  {
    vtkErrorMacro("Invalid dimensions read; check that BytetOrder is set correctly.");
    return -1;
  }

  numPts = dimensions[0] * dimensions[1] * dimensions[2];

  // Skip xCoords
  this->GoldIFile->seekg(sizeof(float) * dimensions[0], ios::cur);
  // Skip yCoords
  this->GoldIFile->seekg(sizeof(float) * dimensions[1], ios::cur);
  // Skip zCoords
  this->GoldIFile->seekg(sizeof(float) * dimensions[2], ios::cur);

  if (iblanked)
  {
    vtkWarningMacro("VTK does not handle blanking for rectilinear grids.");
    this->GoldIFile->seekg(sizeof(int) * numPts, ios::cur);
  }

  // reading next line to check for EOF
  lineRead = this->ReadLine(line);
  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::SkipImageData(char line[256])
{
  char subLine[80];
  int lineRead;
  int iblanked = 0;
  int dimensions[3];
  float origin[3], delta[3];
  int numPts;

  if (sscanf(line, " %*s %*s %s", subLine) == 1)
  {
    if (strncmp(subLine, "iblanked", 8) == 0)
    {
      iblanked = 1;
    }
  }

  this->ReadIntArray(dimensions, 3);
  this->ReadFloatArray(origin, 3);
  this->ReadFloatArray(delta, 3);

  if (iblanked)
  {
    vtkWarningMacro("VTK does not handle blanking for image data.");
    numPts = dimensions[0] * dimensions[1] * dimensions[2];
    if (dimensions[0] < 0 ||
      static_cast<unsigned int>(dimensions[0] * this->SizeOfInt) > this->FileSize ||
      dimensions[1] < 0 ||
      static_cast<unsigned int>(dimensions[1] * this->SizeOfInt) > this->FileSize ||
      dimensions[2] < 0 ||
      static_cast<unsigned int>(dimensions[2] * this->SizeOfInt) > this->FileSize || numPts < 0 ||
      static_cast<unsigned int>(numPts * this->SizeOfInt) > this->FileSize)
    {
      return -1;
    }
    this->GoldIFile->seekg(sizeof(int) * numPts, ios::cur);
  }

  // reading next line to check for EOF
  lineRead = this->ReadLine(line);
  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadMeasuredGeometryFile(
  const char* fileName, int timeStep, vtkMultiBlockDataSet* output)
{
  char line[80], subLine[80];
  vtkIdType i;
  int* pointIds;
  float *xCoords, *yCoords, *zCoords;
  vtkPoints* points = vtkPoints::New();
  vtkPolyData* pd = vtkPolyData::New();

  this->NumberOfNewOutputs++;

  // Initialize
  //
  if (!fileName)
  {
    vtkErrorMacro("A MeasuredFileName must be specified in the case file.");
    return 0;
  }
  std::string sfilename;
  std::string filenameString(fileName);
  this->SanitizeFileName(filenameString);
  if (this->FilePath)
  {
    sfilename = this->FilePath;
    if (sfilename.at(sfilename.length() - 1) != '/')
    {
      sfilename += "/";
    }
    sfilename += filenameString;
    vtkDebugMacro("full path to measured geometry file: " << sfilename);
  }
  else
  {
    sfilename = filenameString;
  }

  if (this->OpenFile(sfilename.c_str()) == 0)
  {
    vtkErrorMacro("Unable to open file: " << sfilename);
    return 0;
  }

  this->ReadLine(line);
  sscanf(line, " %*s %s", subLine);
  if (strncmp(subLine, "Binary", 6) != 0)
  {
    vtkErrorMacro("This is not a binary data set. Try "
      << "vtkEnSightGoldReader.");
    return 0;
  }

  if (this->UseFileSets)
  {
    this->AddFileIndexToCache(fileName);

    i = this->SeekToCachedTimeStep(fileName, timeStep - 1);
    // start w/ the number of TS we skipped, not the one we are at
    // if we are not at the appropriate time step yet, we keep searching
    for (; i < timeStep - 1; i++)
    {
      while (strncmp(line, "BEGIN TIME STEP", 15) != 0)
      {
        this->ReadLine(line);
      }
      // Skip the description line.
      this->ReadLine(line);

      this->ReadLine(line); // "particle coordinates"

      this->ReadInt(&this->NumberOfMeasuredPoints);

      // Skip pointIds
      // this->GoldIFile->ignore(sizeof(int)*this->NumberOfMeasuredPoints);
      // Skip xCoords
      // this->GoldIFile->ignore(sizeof(float)*this->NumberOfMeasuredPoints);
      // Skip yCoords
      // this->GoldIFile->ignore(sizeof(float)*this->NumberOfMeasuredPoints);
      // Skip zCoords
      // this->GoldIFile->ignore(sizeof(float)*this->NumberOfMeasuredPoints);
      this->GoldIFile->seekg(
        (sizeof(float) * 3 + sizeof(int)) * this->NumberOfMeasuredPoints, ios::cur);
      this->ReadLine(line); // END TIME STEP
    }
    while (strncmp(line, "BEGIN TIME STEP", 15) != 0)
    {
      this->ReadLine(line);
    }
    // found a time step -> cache it
    this->AddTimeStepToCache(fileName, i, this->GoldIFile->tellg());
  }

  // Skip the description line.
  this->ReadLine(line);

  this->ReadLine(line); // "particle coordinates"

  this->ReadInt(&this->NumberOfMeasuredPoints);

  pointIds = new int[this->NumberOfMeasuredPoints];
  xCoords = new float[this->NumberOfMeasuredPoints];
  yCoords = new float[this->NumberOfMeasuredPoints];
  zCoords = new float[this->NumberOfMeasuredPoints];
  points->Allocate(this->NumberOfMeasuredPoints);
  pd->AllocateEstimate(this->NumberOfMeasuredPoints, 1);

  // Extract the array of point indices. Note EnSight Manual v8.2 (pp. 559,
  // http://www-vis.lbl.gov/NERSC/Software/ensight/docs82/UserManual.pdf)
  // is wrong in describing the format of binary measured geometry files.
  // As opposed to this description, the actual format employs a 'hybrid'
  // storage scheme. Specifically, point indices are stored in an array,
  // whereas 3D coordinates follow the array in a tuple-by-tuple manner.
  // The following code segment (20+ lines) serves as a fix to bug #9245.
  this->ReadIntArray(pointIds, this->NumberOfMeasuredPoints);

  // Read point coordinates tuple by tuple while each tuple contains three
  // components: (x-cord, y-cord, z-cord)
  int floatSize = sizeof(float);
  for (i = 0; i < this->NumberOfMeasuredPoints; i++)
  {
    this->GoldIFile->read((char*)(xCoords + i), floatSize);
    this->GoldIFile->read((char*)(yCoords + i), floatSize);
    this->GoldIFile->read((char*)(zCoords + i), floatSize);
  }

  if (this->ByteOrder == FILE_LITTLE_ENDIAN)
  {
    vtkByteSwap::Swap4LERange(xCoords, this->NumberOfMeasuredPoints);
    vtkByteSwap::Swap4LERange(yCoords, this->NumberOfMeasuredPoints);
    vtkByteSwap::Swap4LERange(zCoords, this->NumberOfMeasuredPoints);
  }
  else
  {
    vtkByteSwap::Swap4BERange(xCoords, this->NumberOfMeasuredPoints);
    vtkByteSwap::Swap4BERange(yCoords, this->NumberOfMeasuredPoints);
    vtkByteSwap::Swap4BERange(zCoords, this->NumberOfMeasuredPoints);
  }

  // NOTE: EnSight always employs a 1-based indexing scheme and therefore
  // 'if (this->ParticleCoordinatesByIndex)' was removed here. Otherwise
  // the measured geometry could not be proeperly interpreted.
  // This bug was noticed while fixing bug #7453.
  for (i = 0; i < this->NumberOfMeasuredPoints; i++)
  {
    points->InsertNextPoint(xCoords[i], yCoords[i], zCoords[i]);
    pd->InsertNextCell(VTK_VERTEX, 1, &i);
  }

  pd->SetPoints(points);
  this->AddToBlock(output, this->NumberOfGeometryParts, pd);

  points->Delete();
  pd->Delete();
  delete[] pointIds;
  delete[] xCoords;
  delete[] yCoords;
  delete[] zCoords;

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  return 1;
}

//------------------------------------------------------------------------------
bool vtkEnSightGoldBinaryReader::OpenVariableFile(const char* fileName, const char* type)
{
  if (!fileName)
  {
    vtkErrorMacro("Empty " << type << " variable file name.");
    return false;
  }

  std::string sfilename;
  std::string filenameString(fileName);
  this->SanitizeFileName(filenameString);
  if (this->FilePath)
  {
    sfilename = this->FilePath;
    if (sfilename.at(sfilename.length() - 1) != '/')
    {
      sfilename += "/";
    }
    sfilename += filenameString;
    vtkDebugMacro("full path to variable (" << type << ") file: " << sfilename);
  }
  else
  {
    sfilename = filenameString;
  }

  if (this->OpenFile(sfilename.c_str()) == 0)
  {
    vtkErrorMacro("Unable to open file: " << sfilename);
    return false;
  }

  return true;
}

//------------------------------------------------------------------------------
bool vtkEnSightGoldBinaryReader::SkipToTimeStep(const char* fileName, int timeStep,
  vtkMultiBlockDataSet* compositeOutput, int attributeType, int numComponents, bool measured)
{
  if (!this->UseFileSets)
  {
    // nothing to do.
    return true;
  }

  char line[80];

  this->AddFileIndexToCache(fileName);

  int i = this->SeekToCachedTimeStep(fileName, timeStep - 1);
  // start w/ the number of TS we skipped, not the one we are at
  // if we are not at the appropriate time step yet, we keep searching
  for (; i < timeStep - 1; i++)
  {
    this->ReadLine(line);
    while (strncmp(line, "BEGIN TIME STEP", 15) != 0)
    {
      this->ReadLine(line);
    }
    // found a time step -> cache it
    this->AddTimeStepToCache(fileName, i, this->GoldIFile->tellg());

    this->ReadLine(line); // skip the description line

    if (measured)
    {
      auto block = this->GetDataSetFromBlock(compositeOutput, this->NumberOfGeometryParts);
      const auto numElements = block->GetNumberOfElements(attributeType);
      if (numElements)
      {
        this->ReadLine(line);
        // Skip data
        this->GoldIFile->seekg(sizeof(float) * numElements * numComponents, ios::cur);
      }
    }

    int lineRead = this->ReadLine(line); // read "part" line.
    while (lineRead && strncmp(line, "part", 4) == 0)
    {
      int partId;
      this->ReadPartId(&partId);
      partId--; // EnSight starts #ing with 1.
      const auto realId = this->InsertNewPartId(partId);
      auto block = this->GetDataSetFromBlock(compositeOutput, realId);
      const auto numElements = block->GetNumberOfElements(attributeType);
      if (numElements > 0 && attributeType == vtkDataObject::POINT)
      {
        this->ReadLine(line); // "block" or "coordinates"
        // Skip data
        this->GoldIFile->seekg(sizeof(float) * numElements * numComponents, ios::cur);

        // Read the next line.
        lineRead = this->ReadLine(line);
      }
      else if (numElements > 0 && attributeType == vtkDataObject::CELL)
      {
        this->ReadLine(line); // element type or "block"

        // need to find out from CellIds how many cells we have of this
        // element type (and what their ids are) -- IF THIS IS NOT A BLOCK
        // SECTION
        if (strncmp(line, "block", 5) == 0)
        {
          // Skip over float scalars.
          this->GoldIFile->seekg(sizeof(float) * numElements * numComponents, ios::cur);

          // Read the next line.
          lineRead = this->ReadLine(line);
        }
        else
        {
          while (
            lineRead && strncmp(line, "part", 4) != 0 && strncmp(line, "END TIME STEP", 13) != 0)
          {
            const auto elementType = this->GetElementType(line);
            if (elementType == -1)
            {
              vtkErrorMacro("Unknown element type \"" << line << "\"");
              delete this->GoldIFile;
              this->GoldIFile = nullptr;
              return false;
            }
            const auto idx = this->UnstructuredPartIds->IsId(realId);
            const auto numCellsPerElement = this->GetCellIds(idx, elementType)->GetNumberOfIds();
            this->GoldIFile->seekg(sizeof(float) * numCellsPerElement * numComponents, ios::cur);

            // Read the next line.
            lineRead = this->ReadLine(line);
          }
        }
      }
      else // if (numElements == 0)
      {
        // Read the next line.
        lineRead = this->ReadLine(line);
      }
    }
  }

  while (this->ReadLine(line) && strncmp(line, "BEGIN TIME STEP", 15) != 0)
  {
    // keep on advancing;
  }

  return true;
}

//------------------------------------------------------------------------------
bool vtkEnSightGoldBinaryReader::ReadMeasureVariableArray(
  const char* description, vtkMultiBlockDataSet* compositeOutput, int numComponents)
{
  char line[80];
  this->ReadLine(line); // skip the description line

  auto output = static_cast<vtkDataSet*>(
    this->GetDataSetFromBlock(compositeOutput, this->NumberOfGeometryParts));
  const auto numPts = output->GetNumberOfPoints();
  if (numPts)
  {
    // 'this->ReadLine(line)' was removed here, otherwise there would be a
    // problem with timestep retrieval of the measured scalars.
    // This bug was noticed while fixing bug #7453.
    auto array = vtkFloatArray::New();
    array->SetNumberOfComponents(numComponents);
    array->SetNumberOfTuples(numPts);
    this->ReadFloatArray(array->GetPointer(0), numPts * numComponents);
    array->SetName(description);

    auto pd = output->GetPointData();
    pd->AddArray(array);
    if (numComponents == 1 && pd->GetScalars() == nullptr)
    {
      pd->SetScalars(array);
    }
    else if (numComponents == 3 && pd->GetVectors() == nullptr)
    {
      pd->SetVectors(array);
    }
    array->Delete();
  }
  return true;
}

//------------------------------------------------------------------------------
bool vtkEnSightGoldBinaryReader::ReadVariableArray(const char* description,
  vtkMultiBlockDataSet* compositeOutput, int attributeType, int numComponents,
  int component /*=-1*/)
{
  char line[80];

  // read description.
  this->ReadLine(line);

  auto advance = [&line, this]() {
    this->GoldIFile->peek();
    return this->GoldIFile->eof() ? 0 : this->ReadLine(line);
  };

  int lineRead = this->ReadLine(line); // "part".
  while (lineRead && strncmp(line, "part", 4) == 0)
  {
    int partId;
    this->ReadPartId(&partId);
    partId--; // EnSight starts #ing with 1.
    const int realId = this->InsertNewPartId(partId);
    auto output = this->GetDataSetFromBlock(compositeOutput, realId);

    // "part" line can be followed by section identifier "block" or
    // "coordinates" (for nodal data) or "element" (for element data).
    // If the part has zero elements, it may or may not be followed by this section
    // identifier.
    const auto numElements = output->GetNumberOfElements(attributeType);

    lineRead = this->ReadLine(line); // "coordinates", "block", "element" or next "part"
    if (lineRead && strncmp(line, "part", 4) == 0)
    {
      // Part number was not followed by "coordinates" or "block"; we are
      // at the start of another part, skip to next iteration/part to avoid
      // reading anything more.
      continue;
    }

    if (numElements <= 0)
    {
      lineRead = advance();
      continue;
    }

    auto dsa = output->GetAttributes(attributeType);
    vtkSmartPointer<vtkFloatArray> array;

    // For element data (aka cell-data), "part" may be followed by "[element
    // type]" in which case the data is read in chunks rather than a whole.
    if (attributeType != vtkDataObject::CELL || strncmp(line, "block", 5) == 0)
    {
      // read full data.
      array = vtkUtilities::ReadVariableFloats(
        line, this, description, dsa, numElements, numComponents, component);

      lineRead = advance();
    }
    else
    {
      // okay, we are reading in chunks per element type.

      // lets allocate (or get) target array.
      if (component <= 0)
      {
        array = vtk::TakeSmartPointer(vtkFloatArray::New());
        array->SetNumberOfComponents(numComponents);
        array->SetNumberOfTuples(numElements);
      }
      else
      {
        array = vtkFloatArray::SafeDownCast(dsa->GetArray(description));
        assert(array && array->GetNumberOfTuples() == numElements &&
          array->GetNumberOfComponents() == numComponents);
      }

      vtkNew<vtkIdList> srcIds;
      while (lineRead && strncmp(line, "part", 4) != 0 && strncmp(line, "END ", 4) != 0)
      {
        const auto elementType = this->GetElementType(line);
        if (elementType == -1)
        {
          vtkErrorMacro("Unknown element type '" << line << "'");
          return false;
        }

        auto idx = this->UnstructuredPartIds->IsId(realId);
        auto dstIds = this->GetCellIds(idx, elementType);
        const auto numCellsPerElementType = dstIds->GetNumberOfIds();
        auto subarray = vtkUtilities::ReadVariableFloats(
          line, this, description, dsa, numCellsPerElementType, numComponents, component);

        srcIds->SetNumberOfIds(numCellsPerElementType);
        std::iota(srcIds->begin(), srcIds->end(), 0);
        array->InsertTuples(dstIds, srcIds, subarray);

        lineRead = advance();
      }
    }

    if (array)
    {
      array->SetName(description);
      dsa->AddArray(array);
      if (numComponents == 1 && dsa->GetScalars() == nullptr)
      {
        dsa->SetScalars(array);
      }
      else if (numComponents == 3 && dsa->GetVectors() == nullptr)
      {
        dsa->SetVectors(array);
      }
    }
  }

  return true;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadScalarsPerNode(const char* fileName, const char* description,
  int timeStep, vtkMultiBlockDataSet* compositeOutput, int measured, int numberOfComponents,
  int component)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "ScalarsPerNode"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::POINT, 1, measured != 0);

  if (measured)
  {
    this->ReadMeasureVariableArray(description, compositeOutput, 1);
    delete this->GoldIFile;
    this->GoldIFile = nullptr;
    return 1;
  }

  this->ReadVariableArray(
    description, compositeOutput, vtkDataObject::POINT, numberOfComponents, component);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;
  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadAsymmetricTensorsPerNode(const char* fileName,
  const char* description, int timeStep, vtkMultiBlockDataSet* compositeOutput)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "TensorPerNode"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::POINT,
    /* numberOfComponents=*/9, /*measured=*/false);

  this->ReadVariableArray(description, compositeOutput, vtkDataObject::POINT, /*numComponents=*/9);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadVectorsPerNode(const char* fileName, const char* description,
  int timeStep, vtkMultiBlockDataSet* compositeOutput, int measured)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "VectorPerNode"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::POINT,
    /* numberOfComponents=*/3, /*measured=*/measured != 0);

  if (measured)
  {
    this->ReadMeasureVariableArray(description, compositeOutput, 3);
    delete this->GoldIFile;
    this->GoldIFile = nullptr;
    return 1;
  }

  this->ReadVariableArray(description, compositeOutput, vtkDataObject::POINT, 3);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadTensorsPerNode(const char* fileName, const char* description,
  int timeStep, vtkMultiBlockDataSet* compositeOutput)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "TensorsPerNode"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::POINT,
    /* numberOfComponents=*/6, /*measured=*/false);

  this->ReadVariableArray(description, compositeOutput, vtkDataObject::POINT, 6);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadScalarsPerElement(const char* fileName, const char* description,
  int timeStep, vtkMultiBlockDataSet* compositeOutput, int numberOfComponents, int component)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "ScalarPerElement"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::CELL, 1, false);

  this->ReadVariableArray(
    description, compositeOutput, vtkDataObject::CELL, numberOfComponents, component);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;
  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadVectorsPerElement(const char* fileName, const char* description,
  int timeStep, vtkMultiBlockDataSet* compositeOutput)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "VectorPerElement"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::CELL, 3, false);

  this->ReadVariableArray(description, compositeOutput, vtkDataObject::CELL, 3);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;
  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadAsymmetricTensorsPerElement(const char* fileName,
  const char* description, int timeStep, vtkMultiBlockDataSet* compositeOutput)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "AsymmetricTensorPerElement"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::CELL, 9, false);

  this->ReadVariableArray(description, compositeOutput, vtkDataObject::CELL, 9);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::ReadTensorsPerElement(const char* fileName, const char* description,
  int timeStep, vtkMultiBlockDataSet* compositeOutput)
{
  // Initialize
  if (!this->OpenVariableFile(fileName, "TensorPerElement"))
  {
    return 0;
  }

  // If there are multiple timesteps in a file, jump to timestep of interest.
  this->SkipToTimeStep(fileName, timeStep, compositeOutput, vtkDataObject::CELL, 6, false);

  this->ReadVariableArray(description, compositeOutput, vtkDataObject::CELL, 6);

  delete this->GoldIFile;
  this->GoldIFile = nullptr;

  return 1;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::CreateUnstructuredGridOutput(
  int partId, char line[80], const char* name, vtkMultiBlockDataSet* compositeOutput)
{
  int lineRead = 1;
  int i, j;
  int* nodeIdList;
  int numElements;
  int idx, cellId, cellType;
  float *xCoords, *yCoords, *zCoords;

  this->NumberOfNewOutputs++;

  if (this->GetDataSetFromBlock(compositeOutput, partId) == nullptr ||
    !this->GetDataSetFromBlock(compositeOutput, partId)->IsA("vtkUnstructuredGrid"))
  {
    vtkDebugMacro("creating new unstructured output");
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    this->AddToBlock(compositeOutput, partId, ugrid);
    ugrid->Delete();

    this->UnstructuredPartIds->InsertNextId(partId);
  }

  vtkUnstructuredGrid* output =
    vtkUnstructuredGrid::SafeDownCast(this->GetDataSetFromBlock(compositeOutput, partId));
  this->SetBlockName(compositeOutput, partId, name);

  // Clear all cell ids from the last execution, if any.
  idx = this->UnstructuredPartIds->IsId(partId);
  for (i = 0; i < vtkEnSightReader::NUMBER_OF_ELEMENT_TYPES; i++)
  {
    this->GetCellIds(idx, i)->Reset();
  }

  output->Allocate(1000);

  while (lineRead && strncmp(line, "part", 4) != 0)
  {
    if (strncmp(line, "coordinates", 11) == 0)
    {
      vtkDebugMacro("coordinates");
      int numPts;

      this->ReadInt(&numPts);
      if (numPts < 0 || static_cast<unsigned int>(numPts * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro(
          "Invalid number of unstructured points read; check that ByteOrder is set correctly.");
        return -1;
      }

      vtkPoints* points = vtkPoints::New();
      vtkDebugMacro("num. points: " << numPts);

      points->Allocate(numPts);

      if (this->NodeIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numPts + this->FortranSkipBytes, ios::cur);
      }

      xCoords = new float[numPts];
      yCoords = new float[numPts];
      zCoords = new float[numPts];
      this->ReadFloatArray(xCoords, numPts);
      this->ReadFloatArray(yCoords, numPts);
      this->ReadFloatArray(zCoords, numPts);

      for (i = 0; i < numPts; i++)
      {
        points->InsertNextPoint(xCoords[i], yCoords[i], zCoords[i]);
      }

      output->SetPoints(points);
      points->Delete();
      delete[] xCoords;
      delete[] yCoords;
      delete[] zCoords;
    }
    else if (strncmp(line, "point", 5) == 0)
    {
      vtkDebugMacro("point");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of point cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      nodeIdList = new int[numElements];
      this->ReadIntArray(nodeIdList, numElements);

      vtkIdType nodeIds;
      for (i = 0; i < numElements; i++)
      {
        nodeIds = nodeIdList[i] - 1;
        cellId = output->InsertNextCell(VTK_VERTEX, 1, &nodeIds);
        this->GetCellIds(idx, vtkEnSightReader::POINT)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_point", 7) == 0)
    {
      // skipping ghost cells
      vtkDebugMacro("g_point");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of g_point cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      { // skip element ids.
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "bar2", 4) == 0)
    {
      vtkDebugMacro("bar2");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of bar2 cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      nodeIdList = new int[numElements * 2];
      this->ReadIntArray(nodeIdList, numElements * 2);

      vtkIdType nodeIds[2];
      for (i = 0; i < numElements; i++)
      {
        for (j = 0; j < 2; j++)
        {
          nodeIds[j] = nodeIdList[2 * i + j] - 1;
        }
        cellId = output->InsertNextCell(VTK_LINE, 2, nodeIds);
        this->GetCellIds(idx, vtkEnSightReader::BAR2)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_bar2", 6) == 0)
    {
      // skipping ghost cells
      vtkDebugMacro("g_bar2");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of g_bar2 cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * 2 * numElements + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "bar3", 4) == 0)
    {
      vtkDebugMacro("bar3");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of bar3 cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      nodeIdList = new int[numElements * 3];
      this->ReadIntArray(nodeIdList, numElements * 3);

      vtkIdType nodeIds[3];
      for (i = 0; i < numElements; i++)
      {
        nodeIds[0] = nodeIdList[3 * i] - 1;
        nodeIds[1] = nodeIdList[3 * i + 2] - 1;
        nodeIds[2] = nodeIdList[3 * i + 1] - 1;

        cellId = output->InsertNextCell(VTK_QUADRATIC_EDGE, 3, nodeIds);
        this->GetCellIds(idx, vtkEnSightReader::BAR3)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_bar3", 6) == 0)
    {
      // skipping ghost cells
      vtkDebugMacro("g_bar3");

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of g_bar3 cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * 2 * numElements + this->FortranSkipBytes, ios::cur);
    }
    else if (strncmp(line, "nsided", 6) == 0)
    {
      vtkDebugMacro("nsided");
      int* numNodesPerElement;
      int numNodes = 0;
      int nodeCount = 0;

      cellType = vtkEnSightReader::NSIDED;
      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of nsided cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      numNodesPerElement = new int[numElements];
      this->ReadIntArray(numNodesPerElement, numElements);
      for (i = 0; i < numElements; i++)
      {
        numNodes += numNodesPerElement[i];
      }
      nodeIdList = new int[numNodes];
      this->ReadIntArray(nodeIdList, numNodes);

      for (i = 0; i < numElements; i++)
      {
        vtkIdType* nodeIds = new vtkIdType[numNodesPerElement[i]];
        for (j = 0; j < numNodesPerElement[i]; j++)
        {
          nodeIds[j] = nodeIdList[nodeCount] - 1;
          nodeCount++;
        }
        cellId = output->InsertNextCell(VTK_POLYGON, numNodesPerElement[i], nodeIds);
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);

        delete[] nodeIds;
      }

      delete[] nodeIdList;
      delete[] numNodesPerElement;
    }
    else if (strncmp(line, "g_nsided", 8) == 0)
    {
      // skipping ghost cells
      vtkDebugMacro("g_nsided");
      int* numNodesPerElement;
      int numNodes = 0;

      // cellType = vtkEnSightReader::NSIDED;
      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of g_nsided cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      numNodesPerElement = new int[numElements];
      this->ReadIntArray(numNodesPerElement, numElements);
      for (i = 0; i < numElements; i++)
      {
        numNodes += numNodesPerElement[i];
      }
      // Skip nodeIdList.
      this->GoldIFile->seekg(sizeof(int) * numNodes + this->FortranSkipBytes, ios::cur);
      delete[] numNodesPerElement;
    }
    else if (strncmp(line, "tria3", 5) == 0 || strncmp(line, "tria6", 5) == 0)
    {
      if (strncmp(line, "tria6", 5) == 0)
      {
        vtkDebugMacro("tria6");
        cellType = vtkEnSightReader::TRIA6;
      }
      else
      {
        vtkDebugMacro("tria3");
        cellType = vtkEnSightReader::TRIA3;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of triangle cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::TRIA6)
      {
        nodeIdList = new int[numElements * 6];
        this->ReadIntArray(nodeIdList, numElements * 6);
      }
      else
      {
        nodeIdList = new int[numElements * 3];
        this->ReadIntArray(nodeIdList, numElements * 3);
      }

      vtkIdType nodeIds[6];
      for (i = 0; i < numElements; i++)
      {
        if (cellType == vtkEnSightReader::TRIA6)
        {
          for (j = 0; j < 6; j++)
          {
            nodeIds[j] = nodeIdList[6 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUADRATIC_TRIANGLE, 6, nodeIds);
        }
        else
        {
          for (j = 0; j < 3; j++)
          {
            nodeIds[j] = nodeIdList[3 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_TRIANGLE, 3, nodeIds);
        }
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_tria3", 7) == 0 || strncmp(line, "g_tria6", 7) == 0)
    {
      // skipping ghost cells
      if (strncmp(line, "g_tria6", 7) == 0)
      {
        vtkDebugMacro("g_tria6");
        cellType = vtkEnSightReader::TRIA6;
      }
      else
      {
        vtkDebugMacro("g_tria3");
        cellType = vtkEnSightReader::TRIA3;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of triangle cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::TRIA6)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 6 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 3 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "quad4", 5) == 0 || strncmp(line, "quad8", 5) == 0)
    {
      if (strncmp(line, "quad8", 5) == 0)
      {
        vtkDebugMacro("quad8");
        cellType = vtkEnSightReader::QUAD8;
      }
      else
      {
        vtkDebugMacro("quad4");
        cellType = vtkEnSightReader::QUAD4;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of quad cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::QUAD8)
      {
        nodeIdList = new int[numElements * 8];
        this->ReadIntArray(nodeIdList, numElements * 8);
      }
      else
      {
        nodeIdList = new int[numElements * 4];
        this->ReadIntArray(nodeIdList, numElements * 4);
      }

      vtkIdType nodeIds[8];
      for (i = 0; i < numElements; i++)
      {
        if (cellType == vtkEnSightReader::QUAD8)
        {
          for (j = 0; j < 8; j++)
          {
            nodeIds[j] = nodeIdList[8 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUADRATIC_QUAD, 8, nodeIds);
        }
        else
        {
          for (j = 0; j < 4; j++)
          {
            nodeIds[j] = nodeIdList[4 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUAD, 4, nodeIds);
        }
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_quad4", 7) == 0 || strncmp(line, "g_quad8", 7) == 0)
    {
      // skipping ghost cells
      if (strncmp(line, "g_quad8", 7) == 0)
      {
        vtkDebugMacro("g_quad8");
        cellType = vtkEnSightReader::QUAD8;
      }
      else
      {
        vtkDebugMacro("g_quad4");
        cellType = vtkEnSightReader::QUAD4;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of quad cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::QUAD8)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 8 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 4 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }

    else if (strncmp(line, "nfaced", 6) == 0)
    {
      vtkDebugMacro("nfaced");
      int* numFacesPerElement;
      int* numNodesPerFace;
      int* numNodesPerElement;
      int* nodeMarker;
      int numPts = 0;
      int numFaces = 0;
      int numNodes = 0;
      int faceCount = 0;
      int nodeCount = 0;
      int elementNodeCount = 0;

      cellType = vtkEnSightReader::NFACED;
      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of nfaced cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      // array: number of faces per element
      numFacesPerElement = new int[numElements];
      this->ReadIntArray(numFacesPerElement, numElements);

      // array: number of nodes per face
      for (i = 0; i < numElements; i++)
      {
        numFaces += numFacesPerElement[i];
      }
      numNodesPerFace = new int[numFaces];
      this->ReadIntArray(numNodesPerFace, numFaces);

      // array: number of nodes per element
      // number of faces of all elements
      numNodesPerElement = new int[numElements];
      for (i = 0; i < numElements; i++)
      {
        numNodesPerElement[i] = 0;
        for (j = 0; j < numFacesPerElement[i]; j++)
        {
          numNodesPerElement[i] += numNodesPerFace[faceCount + j];
        }
        faceCount += numFacesPerElement[i];
      }

      // xxx begin
      // delete [] numFacesPerElement;
      // delete [] numNodesPerFace;
      // xxx end

      // number of nodes of all elements
      for (i = 0; i < numElements; i++)
      {
        numNodes += numNodesPerElement[i];
      }

      // allocate and init markers to determine unique points
      numPts = output->GetNumberOfPoints();
      nodeMarker = new int[numPts];
      for (i = 0; i < numPts; i++)
      {
        nodeMarker[i] = -1;
      }

      // array: node Ids of all elements
      // NOTE:  each node Id is usually referenced multiple times in a
      //        polyhedron and therefore nodeIdList is not a set of
      //        UNIQUE point Ids (instead it an RAW list)
      nodeIdList = new int[numNodes];
      this->ReadIntArray(nodeIdList, numNodes);

      // yyy begin
      int k;                        // indexing each node Id of a face
      int faceIdx = 0;              // indexing faces throughout all polyhedra
      int nodeIdx = 0;              // indexing nodes throughout all polyhedra
      int arayIdx = 0;              // indexing the array of Ids (info of faces)
      vtkIdType* faceAry = nullptr; // array of Ids describing a vtkPolyhedron
      // yyy end

      for (i = 0; i < numElements; i++)
      {
        elementNodeCount = 0;
        vtkIdType* nodeIds = new vtkIdType[numNodesPerElement[i]];

        // yyy begin
        arayIdx = 0;
        faceAry = new vtkIdType[numFacesPerElement[i] + numNodesPerElement[i]];
        for (j = 0; j < numFacesPerElement[i]; j++, faceIdx++)
        {
          faceAry[arayIdx++] = numNodesPerFace[faceIdx];

          for (k = 0; k < numNodesPerFace[faceIdx]; k++)
          {
            faceAry[arayIdx++] = nodeIdList[nodeIdx++] - 1;
          }
        }
        // yyy end

        for (j = 0; j < numNodesPerElement[i]; j++)
        {
          if (nodeMarker[nodeIdList[nodeCount] - 1] < i)
          {
            nodeIds[elementNodeCount] = nodeIdList[nodeCount] - 1;
            nodeMarker[nodeIdList[nodeCount] - 1] = i;
            elementNodeCount += 1;
          }
          nodeCount++;
        }

        // xxx begin
        // cellId = output->InsertNextCell( VTK_CONVEX_POINT_SET,
        //                                 elementNodeCount, nodeIds );
        // xxx end

        // yyy begin
        cellId = output->InsertNextCell(
          VTK_POLYHEDRON, elementNodeCount, nodeIds, numFacesPerElement[i], faceAry);
        delete[] faceAry;
        faceAry = nullptr;
        // yyy end

        this->GetCellIds(idx, cellType)->InsertNextId(cellId);

        delete[] nodeIds;
      }

      // yyy begin
      delete[] numNodesPerFace;
      delete[] numFacesPerElement;
      numNodesPerFace = nullptr;
      numFacesPerElement = nullptr;
      // yyy end

      delete[] nodeMarker;
      delete[] nodeIdList;
      delete[] numNodesPerElement;
    }
    else if (strncmp(line, "tetra4", 6) == 0 || strncmp(line, "tetra10", 7) == 0)
    {
      if (strncmp(line, "tetra10", 7) == 0)
      {
        vtkDebugMacro("tetra10");
        cellType = vtkEnSightReader::TETRA10;
      }
      else
      {
        vtkDebugMacro("tetra4");
        cellType = vtkEnSightReader::TETRA4;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro(
          "Invalid number of tetrahedral cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::TETRA10)
      {
        nodeIdList = new int[numElements * 10];
        this->ReadIntArray(nodeIdList, numElements * 10);
      }
      else
      {
        nodeIdList = new int[numElements * 4];
        this->ReadIntArray(nodeIdList, numElements * 4);
      }

      vtkIdType nodeIds[10];
      for (i = 0; i < numElements; i++)
      {
        if (cellType == vtkEnSightReader::TETRA10)
        {
          for (j = 0; j < 10; j++)
          {
            nodeIds[j] = nodeIdList[10 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUADRATIC_TETRA, 10, nodeIds);
        }
        else
        {
          for (j = 0; j < 4; j++)
          {
            nodeIds[j] = nodeIdList[4 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_TETRA, 4, nodeIds);
        }
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_tetra4", 8) == 0 || strncmp(line, "g_tetra10", 9) == 0)
    {
      // skipping ghost cells
      if (strncmp(line, "g_tetra10", 9) == 0)
      {
        vtkDebugMacro("g_tetra10");
        cellType = vtkEnSightReader::TETRA10;
      }
      else
      {
        vtkDebugMacro("g_tetra4");
        cellType = vtkEnSightReader::TETRA4;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro(
          "Invalid number of tetrahedral cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::TETRA10)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 10 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 4 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "pyramid5", 8) == 0 || strncmp(line, "pyramid13", 9) == 0)
    {
      if (strncmp(line, "pyramid13", 9) == 0)
      {
        vtkDebugMacro("pyramid13");
        cellType = vtkEnSightReader::PYRAMID13;
      }
      else
      {
        vtkDebugMacro("pyramid5");
        cellType = vtkEnSightReader::PYRAMID5;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of pyramid cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::PYRAMID13)
      {
        nodeIdList = new int[numElements * 13];
        this->ReadIntArray(nodeIdList, numElements * 13);
      }
      else
      {
        nodeIdList = new int[numElements * 5];
        this->ReadIntArray(nodeIdList, numElements * 5);
      }

      vtkIdType nodeIds[13];
      for (i = 0; i < numElements; i++)
      {
        if (cellType == vtkEnSightReader::PYRAMID13)
        {
          for (j = 0; j < 13; j++)
          {
            nodeIds[j] = nodeIdList[13 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUADRATIC_PYRAMID, 13, nodeIds);
        }
        else
        {
          for (j = 0; j < 5; j++)
          {
            nodeIds[j] = nodeIdList[5 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_PYRAMID, 5, nodeIds);
        }
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_pyramid5", 10) == 0 || strncmp(line, "g_pyramid13", 11) == 0)
    {
      // skipping ghost cells
      if (strncmp(line, "g_pyramid13", 11) == 0)
      {
        vtkDebugMacro("g_pyramid13");
        cellType = vtkEnSightReader::PYRAMID13;
      }
      else
      {
        vtkDebugMacro("g_pyramid5");
        cellType = vtkEnSightReader::PYRAMID5;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of pyramid cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::PYRAMID13)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 13 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 5 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "hexa8", 5) == 0 || strncmp(line, "hexa20", 6) == 0)
    {
      if (strncmp(line, "hexa20", 6) == 0)
      {
        vtkDebugMacro("hexa20");
        cellType = vtkEnSightReader::HEXA20;
      }
      else
      {
        vtkDebugMacro("hexa8");
        cellType = vtkEnSightReader::HEXA8;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of hexahedral cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::HEXA20)
      {
        nodeIdList = new int[numElements * 20];
        this->ReadIntArray(nodeIdList, numElements * 20);
      }
      else
      {
        nodeIdList = new int[numElements * 8];
        this->ReadIntArray(nodeIdList, numElements * 8);
      }

      vtkIdType nodeIds[20];
      for (i = 0; i < numElements; i++)
      {
        if (cellType == vtkEnSightReader::HEXA20)
        {
          for (j = 0; j < 20; j++)
          {
            nodeIds[j] = nodeIdList[20 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUADRATIC_HEXAHEDRON, 20, nodeIds);
        }
        else
        {
          for (j = 0; j < 8; j++)
          {
            nodeIds[j] = nodeIdList[8 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_HEXAHEDRON, 8, nodeIds);
        }
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_hexa8", 7) == 0 || strncmp(line, "g_hexa20", 8) == 0)
    {
      // skipping ghost cells
      if (strncmp(line, "g_hexa20", 8) == 0)
      {
        vtkDebugMacro("g_hexa20");
        cellType = vtkEnSightReader::HEXA20;
      }
      else
      {
        vtkDebugMacro("g_hexa8");
        cellType = vtkEnSightReader::HEXA8;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of hexahedral cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::HEXA20)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 20 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 8 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "penta6", 6) == 0 || strncmp(line, "penta15", 7) == 0)
    {
      if (strncmp(line, "penta15", 7) == 0)
      {
        vtkDebugMacro("penta15");
        cellType = vtkEnSightReader::PENTA15;
      }
      else
      {
        vtkDebugMacro("penta6");
        cellType = vtkEnSightReader::PENTA6;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of pentagonal cells; check that ByteOrder is set correctly.");
        return -1;
      }

      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::PENTA15)
      {
        nodeIdList = new int[numElements * 15];
        this->ReadIntArray(nodeIdList, numElements * 15);
      }
      else
      {
        nodeIdList = new int[numElements * 6];
        this->ReadIntArray(nodeIdList, numElements * 6);
      }

      const unsigned char penta6Map[6] = { 0, 2, 1, 3, 5, 4 };
      const unsigned char penta15Map[15] = { 0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13 };

      vtkIdType nodeIds[15];
      for (i = 0; i < numElements; i++)
      {
        if (cellType == vtkEnSightReader::PENTA15)
        {
          for (j = 0; j < 15; j++)
          {
            nodeIds[penta15Map[j]] = nodeIdList[15 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_QUADRATIC_WEDGE, 15, nodeIds);
        }
        else
        {
          for (j = 0; j < 6; j++)
          {
            nodeIds[penta6Map[j]] = nodeIdList[6 * i + j] - 1;
          }
          cellId = output->InsertNextCell(VTK_WEDGE, 6, nodeIds);
        }
        this->GetCellIds(idx, cellType)->InsertNextId(cellId);
      }

      delete[] nodeIdList;
    }
    else if (strncmp(line, "g_penta6", 8) == 0 || strncmp(line, "g_penta15", 9) == 0)
    {
      // skipping ghost cells
      if (strncmp(line, "g_penta15", 9) == 0)
      {
        vtkDebugMacro("g_penta15");
        cellType = vtkEnSightReader::PENTA15;
      }
      else
      {
        vtkDebugMacro("g_penta6");
        cellType = vtkEnSightReader::PENTA6;
      }

      this->ReadInt(&numElements);
      if (numElements < 0 ||
        static_cast<unsigned int>(numElements * this->SizeOfInt) > this->FileSize)
      {
        vtkErrorMacro("Invalid number of pentagonal cells; check that ByteOrder is set correctly.");
        return -1;
      }
      if (this->ElementIdsListed)
      {
        this->GoldIFile->seekg(sizeof(int) * numElements + this->FortranSkipBytes, ios::cur);
      }

      if (cellType == vtkEnSightReader::PENTA15)
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 15 * numElements + this->FortranSkipBytes, ios::cur);
      }
      else
      {
        // Skip nodeIdList.
        this->GoldIFile->seekg(sizeof(int) * 6 * numElements + this->FortranSkipBytes, ios::cur);
      }
    }
    else if (strncmp(line, "END TIME STEP", 13) == 0)
    {
      return 1;
    }
    else if (this->IS && this->IS->fail())
    {
      // May want consistency check here?
      // vtkWarningMacro("EOF on geometry file");
      return 1;
    }
    else
    {
      vtkErrorMacro("undefined geometry file line '" << line << "'");
      return -1;
    }
    this->GoldIFile->peek();
    if (this->GoldIFile->eof())
    {
      lineRead = 0;
      continue;
    }
    lineRead = this->ReadLine(line);
  }
  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::CreateStructuredGridOutput(
  int partId, char line[80], const char* name, vtkMultiBlockDataSet* compositeOutput)
{
  char subLine[80];
  int lineRead;
  int iblanked = 0;
  int dimensions[3];
  int i;
  vtkPoints* points = vtkPoints::New();
  int numPts;
  float *xCoords, *yCoords, *zCoords;

  this->NumberOfNewOutputs++;

  vtkDataSet* ds = this->GetDataSetFromBlock(compositeOutput, partId);
  if (ds == nullptr || !ds->IsA("vtkStructuredGrid"))
  {
    vtkDebugMacro("creating new structured grid output");
    vtkStructuredGrid* sgrid = vtkStructuredGrid::New();
    this->AddToBlock(compositeOutput, partId, sgrid);
    sgrid->Delete();
    ds = sgrid;
  }

  vtkStructuredGrid* output = vtkStructuredGrid::SafeDownCast(ds);
  this->SetBlockName(compositeOutput, partId, name);

  if (sscanf(line, " %*s %s", subLine) == 1)
  {
    if (strncmp(subLine, "iblanked", 8) == 0)
    {
      iblanked = 1;
    }
  }

  this->ReadIntArray(dimensions, 3);
  numPts = dimensions[0] * dimensions[1] * dimensions[2];
  if (dimensions[0] < 0 ||
    static_cast<unsigned int>(dimensions[0] * this->SizeOfInt) > this->FileSize ||
    dimensions[1] < 0 ||
    static_cast<unsigned int>(dimensions[1] * this->SizeOfInt) > this->FileSize ||
    dimensions[2] < 0 ||
    static_cast<unsigned int>(dimensions[2] * this->SizeOfInt) > this->FileSize || numPts < 0 ||
    static_cast<unsigned int>(numPts * this->SizeOfInt) > this->FileSize)
  {
    vtkErrorMacro("Invalid dimensions read; check that ByteOrder is set correctly.");
    points->Delete();
    return -1;
  }
  output->SetDimensions(dimensions);
  points->Allocate(numPts);

  xCoords = new float[numPts];
  yCoords = new float[numPts];
  zCoords = new float[numPts];
  this->ReadFloatArray(xCoords, numPts);
  this->ReadFloatArray(yCoords, numPts);
  this->ReadFloatArray(zCoords, numPts);

  for (i = 0; i < numPts; i++)
  {
    points->InsertNextPoint(xCoords[i], yCoords[i], zCoords[i]);
  }
  output->SetPoints(points);
  if (iblanked)
  {
    int* iblanks = new int[numPts];
    this->ReadIntArray(iblanks, numPts);

    for (i = 0; i < numPts; i++)
    {
      if (!iblanks[i])
      {
        output->BlankPoint(i);
      }
    }
    delete[] iblanks;
  }

  points->Delete();
  delete[] xCoords;
  delete[] yCoords;
  delete[] zCoords;

  this->GoldIFile->peek();
  if (this->GoldIFile->eof())
  {
    lineRead = 0;
  }
  else
  {
    lineRead = this->ReadLine(line);
  }

  if (strncmp(line, "node_ids", 8) == 0)
  {
    int* nodeIds = new int[numPts];
    this->ReadIntArray(nodeIds, numPts);
    lineRead = this->ReadLine(line);
    delete[] nodeIds;
  }
  if (strncmp(line, "element_ids", 11) == 0)
  {
    int numElements = (dimensions[0] - 1) * (dimensions[1] - 1) * (dimensions[2] - 1);
    int* elementIds = new int[numElements];
    this->ReadIntArray(elementIds, numElements);
    lineRead = this->ReadLine(line);
    delete[] elementIds;
  }

  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::CreateRectilinearGridOutput(
  int partId, char line[80], const char* name, vtkMultiBlockDataSet* compositeOutput)
{
  char subLine[80];
  int lineRead;
  int iblanked = 0;
  int dimensions[3];
  int i;
  vtkFloatArray* xCoords = vtkFloatArray::New();
  vtkFloatArray* yCoords = vtkFloatArray::New();
  vtkFloatArray* zCoords = vtkFloatArray::New();
  float* tempCoords;
  int numPts;

  this->NumberOfNewOutputs++;

  vtkDataSet* ds = this->GetDataSetFromBlock(compositeOutput, partId);
  if (ds == nullptr || !ds->IsA("vtkRectilinearGrid"))
  {
    vtkDebugMacro("creating new rectilinear grid output");
    vtkRectilinearGrid* rgrid = vtkRectilinearGrid::New();
    this->AddToBlock(compositeOutput, partId, rgrid);
    rgrid->Delete();
    ds = rgrid;
  }

  vtkRectilinearGrid* output = vtkRectilinearGrid::SafeDownCast(ds);

  this->SetBlockName(compositeOutput, partId, name);

  if (sscanf(line, " %*s %*s %s", subLine) == 1)
  {
    if (strncmp(subLine, "iblanked", 8) == 0)
    {
      iblanked = 1;
    }
  }

  this->ReadIntArray(dimensions, 3);
  if (dimensions[0] < 0 ||
    static_cast<unsigned int>(dimensions[0] * this->SizeOfInt) > this->FileSize ||
    dimensions[1] < 0 ||
    static_cast<unsigned int>(dimensions[1] * this->SizeOfInt) > this->FileSize ||
    dimensions[2] < 0 ||
    static_cast<unsigned int>(dimensions[2] * this->SizeOfInt) > this->FileSize ||
    (dimensions[0] + dimensions[1] + dimensions[2]) < 0 ||
    static_cast<unsigned int>((dimensions[0] + dimensions[1] + dimensions[2]) * this->SizeOfInt) >
      this->FileSize)
  {
    vtkErrorMacro("Invalid dimensions read; check that BytetOrder is set correctly.");
    xCoords->Delete();
    yCoords->Delete();
    zCoords->Delete();
    return -1;
  }

  output->SetDimensions(dimensions);
  xCoords->Allocate(dimensions[0]);
  yCoords->Allocate(dimensions[1]);
  zCoords->Allocate(dimensions[2]);

  tempCoords = new float[dimensions[0]];
  this->ReadFloatArray(tempCoords, dimensions[0]);
  for (i = 0; i < dimensions[0]; i++)
  {
    xCoords->InsertNextTuple(&tempCoords[i]);
  }
  delete[] tempCoords;
  tempCoords = new float[dimensions[1]];
  this->ReadFloatArray(tempCoords, dimensions[1]);
  for (i = 0; i < dimensions[1]; i++)
  {
    yCoords->InsertNextTuple(&tempCoords[i]);
  }
  delete[] tempCoords;
  tempCoords = new float[dimensions[2]];
  this->ReadFloatArray(tempCoords, dimensions[2]);
  for (i = 0; i < dimensions[2]; i++)
  {
    zCoords->InsertNextTuple(&tempCoords[i]);
  }
  delete[] tempCoords;
  if (iblanked)
  {
    vtkWarningMacro("VTK does not handle blanking for rectilinear grids.");
    numPts = dimensions[0] * dimensions[1] * dimensions[2];
    int* tempArray = new int[numPts];
    this->ReadIntArray(tempArray, numPts);
    delete[] tempArray;
  }

  output->SetXCoordinates(xCoords);
  output->SetYCoordinates(yCoords);
  output->SetZCoordinates(zCoords);

  xCoords->Delete();
  yCoords->Delete();
  zCoords->Delete();

  // reading next line to check for EOF
  lineRead = this->ReadLine(line);
  return lineRead;
}

//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::CreateImageDataOutput(
  int partId, char line[80], const char* name, vtkMultiBlockDataSet* compositeOutput)
{
  char subLine[80];
  int lineRead;
  int iblanked = 0;
  int dimensions[3];
  float origin[3], delta[3];
  int numPts;

  this->NumberOfNewOutputs++;

  vtkDataSet* ds = this->GetDataSetFromBlock(compositeOutput, partId);
  if (ds == nullptr || !ds->IsA("vtkImageData"))
  {
    vtkDebugMacro("creating new image data output");
    vtkImageData* idata = vtkImageData::New();
    this->AddToBlock(compositeOutput, partId, idata);
    idata->Delete();
    ds = idata;
  }

  vtkImageData* output = vtkImageData::SafeDownCast(ds);

  this->SetBlockName(compositeOutput, partId, name);

  if (sscanf(line, " %*s %*s %s", subLine) == 1)
  {
    if (strncmp(subLine, "iblanked", 8) == 0)
    {
      iblanked = 1;
    }
  }

  this->ReadIntArray(dimensions, 3);
  output->SetDimensions(dimensions);
  this->ReadFloatArray(origin, 3);
  output->SetOrigin(origin[0], origin[1], origin[2]);
  this->ReadFloatArray(delta, 3);
  output->SetSpacing(delta[0], delta[1], delta[2]);

  if (iblanked)
  {
    vtkWarningMacro("VTK does not handle blanking for image data.");
    numPts = dimensions[0] * dimensions[1] * dimensions[2];
    if (dimensions[0] < 0 ||
      static_cast<unsigned int>(dimensions[0] * this->SizeOfInt) > this->FileSize ||
      dimensions[1] < 0 ||
      static_cast<unsigned int>(dimensions[1] * this->SizeOfInt) > this->FileSize ||
      dimensions[2] < 0 ||
      static_cast<unsigned int>(dimensions[2] * this->SizeOfInt) > this->FileSize || numPts < 0 ||
      static_cast<unsigned int>(numPts * this->SizeOfInt) > this->FileSize)
    {
      return -1;
    }
    int* tempArray = new int[numPts];
    this->ReadIntArray(tempArray, numPts);
    delete[] tempArray;
  }

  // reading next line to check for EOF
  lineRead = this->ReadLine(line);
  return lineRead;
}

// Internal function to read in a line up to 80 characters.
// Returns zero if there was an error.
int vtkEnSightGoldBinaryReader::ReadLine(char result[80])
{
  if (!this->GoldIFile->read(result, 80))
  {
    // The read fails when reading the last part/array when there are no points.
    // I took out the error macro as a temporary fix.
    // We need to determine what EnSight does when the part with zero point
    // is not the last, and change the read array method.
    // int fixme; // I do not have a file to test with yet.
    vtkDebugMacro("Read failed");
    return 0;
  }
  // fix to the memory leakage problem detected by Valgrind
  result[79] = '\0';

  if (this->Fortran)
  {
    // strncpy cannot be used for overlapping buffers
    int i = 0;
    for (; i < 76 && result[i + 4] != '\0'; ++i)
    {
      result[i] = result[i + 4];
    }
    for (; i < 76; ++i)
    {
      result[i] = '\0';
    }
    result[76] = 0;
    // better read an extra 8 bytes to prevent error next time
    char dummy[8];
    if (!this->GoldIFile->read(dummy, 8))
    {
      vtkDebugMacro("Read (fortran) failed");
      return 0;
    }
  }

  return 1;
}

// Internal function to read a single integer.
// Returns zero if there was an error.
// Sets byte order so that part id is reasonable.
int vtkEnSightGoldBinaryReader::ReadPartId(int* result)
{
  // first swap like normal.
  if (this->ReadInt(result) == 0)
  {
    vtkErrorMacro("Read failed");
    return 0;
  }

  // second: try an experimental byte swap.
  // Only experiment if byte order is not set.
  if (this->ByteOrder == FILE_UNKNOWN_ENDIAN)
  {
    int tmpLE = *result;
    int tmpBE = *result;
    vtkByteSwap::Swap4LE(&tmpLE);
    vtkByteSwap::Swap4BE(&tmpBE);

    if (tmpLE >= 0 && tmpLE < MAXIMUM_PART_ID)
    {
      this->ByteOrder = FILE_LITTLE_ENDIAN;
      *result = tmpLE;
      return 1;
    }
    if (tmpBE >= 0 && tmpBE < MAXIMUM_PART_ID)
    {
      this->ByteOrder = FILE_BIG_ENDIAN;
      *result = tmpBE;
      return 1;
    }
    vtkErrorMacro("Byte order could not be determined.");
    return 0;
  }

  return 1;
}

// Internal function to read a single integer.
// Returns zero if there was an error.
int vtkEnSightGoldBinaryReader::ReadInt(int* result)
{
  char dummy[4];
  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  if (!this->GoldIFile->read((char*)result, sizeof(int)))
  {
    vtkErrorMacro("Read failed");
    return 0;
  }

  if (this->ByteOrder == FILE_LITTLE_ENDIAN)
  {
    vtkByteSwap::Swap4LE(result);
  }
  else if (this->ByteOrder == FILE_BIG_ENDIAN)
  {
    vtkByteSwap::Swap4BE(result);
  }

  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  return 1;
}

// Internal function to read a single float.
// Returns zero if there was an error.
int vtkEnSightGoldBinaryReader::ReadFloat(float* result)
{
  char dummy[4];
  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  if (!this->GoldIFile->read((char*)result, sizeof(float)))
  {
    vtkErrorMacro("Read failed");
    return 0;
  }

  if (this->ByteOrder == FILE_LITTLE_ENDIAN)
  {
    vtkByteSwap::Swap4LE(result);
  }
  else if (this->ByteOrder == FILE_BIG_ENDIAN)
  {
    vtkByteSwap::Swap4BE(result);
  }

  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  return 1;
}

// Internal function to read an integer array.
// Returns zero if there was an error.
int vtkEnSightGoldBinaryReader::ReadIntArray(int* result, int numInts)
{
  if (numInts <= 0)
  {
    return 1;
  }

  char dummy[4];
  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  if (!this->GoldIFile->read((char*)result, sizeof(int) * numInts))
  {
    vtkErrorMacro("Read failed.");
    return 0;
  }

  if (this->ByteOrder == FILE_LITTLE_ENDIAN)
  {
    vtkByteSwap::Swap4LERange(result, numInts);
  }
  else
  {
    vtkByteSwap::Swap4BERange(result, numInts);
  }

  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  return 1;
}

// Internal function to read a single vtkTypeInt64.
// Returns zero if there was an error.
int vtkEnSightGoldBinaryReader::ReadLong(vtkTypeInt64* result)
{
  char dummy[4];
  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  if (!this->GoldIFile->read((char*)result, sizeof(vtkTypeInt64)))
  {
    vtkErrorMacro("Read failed");
    return 0;
  }

  if (this->ByteOrder == FILE_LITTLE_ENDIAN)
  {
    vtkByteSwap::Swap8LE(result);
  }
  else if (this->ByteOrder == FILE_BIG_ENDIAN)
  {
    vtkByteSwap::Swap8BE(result);
  }

  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  return 1;
}

// Internal function to read a float array.
// Returns zero if there was an error.
int vtkEnSightGoldBinaryReader::ReadFloatArray(float* result, int numFloats)
{
  if (numFloats <= 0)
  {
    return 1;
  }

  char dummy[4];
  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }

  if (!this->GoldIFile->read((char*)result, sizeof(float) * numFloats))
  {
    vtkErrorMacro("Read failed");
    return 0;
  }

  if (this->ByteOrder == FILE_LITTLE_ENDIAN)
  {
    vtkByteSwap::Swap4LERange(result, numFloats);
  }
  else
  {
    vtkByteSwap::Swap4BERange(result, numFloats);
  }

  if (this->Fortran)
  {
    if (!this->GoldIFile->read(dummy, 4))
    {
      vtkErrorMacro("Read (fortran) failed.");
      return 0;
    }
  }
  return 1;
}

//------------------------------------------------------------------------------
void vtkEnSightGoldBinaryReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// Seeks the IFile to the cached timestep nearest the target timestep.
// Returns the actually sought to timestep
//------------------------------------------------------------------------------
int vtkEnSightGoldBinaryReader::SeekToCachedTimeStep(const char* fileName, int realTimeStep)
{
  typedef vtkEnSightGoldBinaryReader::FileOffsetMapInternal MapType;
  typedef MapType::const_iterator MapNameIterator;
  typedef MapType::value_type::second_type::const_iterator FileOffsetIterator;

  int j = 0;

  MapNameIterator nameIterator = this->FileOffsets->Map.find(fileName);
  if (nameIterator == this->FileOffsets->Map.end())
  {
    return j;
  }

  // Try to find the nearest time step for which we know the offset
  for (int i = realTimeStep; i >= 0; i--)
  {
    FileOffsetIterator fileOffsetIterator = nameIterator->second.find(i);
    if (fileOffsetIterator != nameIterator->second.end())
    {
      // we need to account for the last 80 characters as where we need to seek,
      // as we need to be at the BEGIN TIMESTEP keyword and not
      // the description line
      this->GoldIFile->seekg(fileOffsetIterator->second - 80l, ios::beg);
      j = i;
      break;
    }
  }
  return j;
}

//------------------------------------------------------------------------------
// Add a cached file offset
void vtkEnSightGoldBinaryReader::AddTimeStepToCache(
  const char* fileName, int realTimeStep, vtkTypeInt64 address)
{
  if (this->FileOffsets->Map.find(fileName) == this->FileOffsets->Map.end())
  {
    std::map<int, vtkTypeInt64> tsMap;
    this->FileOffsets->Map[fileName] = tsMap;
  }
  this->FileOffsets->Map[fileName][realTimeStep] = address;
}

//------------------------------------------------------------------------------
void vtkEnSightGoldBinaryReader::AddFileIndexToCache(const char* fileName)
{
  // only read the file index if we have not searched for the file index before
  if (this->FileOffsets->Map.find(fileName) == this->FileOffsets->Map.end())
  {
    char line[80];
    vtkTypeInt64 addr;
    int numTS;

    // We add an empty map to prevent further attempts at reading the file index
    std::map<int, vtkTypeInt64> tsMap;
    this->FileOffsets->Map[fileName] = tsMap;

    // Read the last 80 characters (+ a vtkTypeInt64) of the file and check for FILE_INDEX
    vtkIdType seekOffset = (vtkIdType(-80) * static_cast<vtkIdType>(sizeof(char))) -
      static_cast<vtkIdType>(sizeof(vtkTypeInt64));
    this->GoldIFile->seekg(seekOffset, ios::end);

    // right before the FILE_INDEX entry we might find the address of the index start
    this->ReadLong(&addr);
    this->ReadLine(line);

    if (strncmp(line, "FILE_INDEX", 10) == 0)
    {
      // jump to beginning of the index and add all time steps
      this->GoldIFile->seekg(addr, ios::beg);
      this->ReadInt(&numTS);

      for (int i = 0; i < numTS; ++i)
      {
        this->ReadLong(&addr);
        // The file index points at the description line, while VTK points at BEGIN TIMESTEP
        this->FileOffsets->Map[fileName][i] = addr;
      }
    }
  }
  this->GoldIFile->seekg(0l, ios::beg);
}
VTK_ABI_NAMESPACE_END

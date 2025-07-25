
## Test-dummyLib  (directory: 00-dummy)

    Minimal compilation test with wmake, without OpenFOAM libraries.

    The application and library can also serve as a minimal test case for
    wmake, or to provide a minimal library/executable target for testing.


## Test-machine-sizes  (directory: 00-machine-sizes)

    Test the sizeof for basic types.
    Can be compiled and run without any OpenFOAM libraries.

        g++ -std=c++11 -oTest-machine-sizes Test-machine-sizes.cpp


## Test-openmp  (directory: 00-openmp)

    Simple test program for compiling/running openmp


## Test-BinSum  (directory: BinSum)

    Test BinSum container


## Test-CircularBuffer  (directory: CircularBuffer)

    Basic tests for CircularBuffer behaviour and characteristics


## Test-Circulator  (directory: Circulator)

    Tests for Circulator and ConstCirculator


## Test-CompactIOList  (directory: CompactIOList)

    Simple demonstration and test application for the CompactIOList container


## Test-CompactListList  (directory: CompactListList)

    Simple demonstration and test application for the CompactListList class.


## Test-DLList  (directory: DLList)

    Tests for doubly-linked lists


## Test-DiagTensor  (directory: DiagTensor)

    Tests for \c DiagTensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-Dictionary  (directory: Dictionary)

    Tests for Dictionary (not dictionary)


## Test-DirLister  (directory: DirLister)

    Test functionality of DirLister


## Test-Distribution  (directory: Distribution)

    Test the Distribution class

    Plot normal distribution test in gnuplot using:

    \verbatim
    normalDistribution(mean, sigma, x) = \
        sqrt(1.0/(2.0*pi*sigma**2))*exp(-(x - mean)**2.0/(2.0*sigma**2))

    plot normalDistribution(8.5, 2.5, x), "Distribution_scalar_test_x" w p
    \endverbatim


## Test-DynamicField  (directory: DynamicField)

    Tests for DynamicField base functionality


## Test-DynamicList  (directory: DynamicList)

    Tests for DynamicList base functionality


## Test-DynamicList2  (directory: DynamicList2)

    Test allocation patterns when reading into an existing list.


## Test-Enum  (directory: Enum)

    Testing of Enum lookups.


## Test-FieldFields1  (directory: FieldFields1)

- no description

## Test-FieldFields2  (directory: FieldFields2)

- no description

## Test-FixedList  (directory: FixedList)

    Simple tests and examples for FixedList

See also
    Foam::FixedList


## Test-FixedList2  (directory: FixedList2)

    Test speeds, usability of some List/FixedList operations

See also
    Foam::FixedList


## Test-Function1  (directory: Function1)

    Tests Function1


## Test-GAMGAgglomeration  (directory: GAMGAgglomeration)

    Test application for GAMG agglomeration. Hardcoded to expect GAMG on p.


## Test-HashPtrTable  (directory: HashPtrTable)

    Tests for HashPtrTable base functionality


## Test-hashSet  (directory: HashSet)

    Some simple HashSet tests


## Test-HashTable1  (directory: HashTable1)

- no description

## Test-HashTable2  (directory: HashTable2)

    Miscellaneous tests for HashTable


## Test-HashTable3  (directory: HashTable3)

    Test speeds for some HashTable operations


## Test-HashTable4  (directory: HashTable4)

    Test HashTable resizing


## Test-Hashing1  (directory: Hashing1)

    Test/verify overloads of Hash function


## Test-Hashing2  (directory: Hashing2)



## Test-HashingSpeed  (directory: HashingSpeed)

- no description

## Test-ICharStream1  (directory: ICharStream1)



## Test-IFstream  (directory: IFstream)

    Tests on low-level reading


## Test-IOField  (directory: IOField)

    Test the processor-local reading of IOField (used in the lagrangian libs)


## Test-IOobjectList  (directory: IOobjectList)

    Basic tests of IOobjectList


## Test-ISLList  (directory: ISLList)



## Test-IStringStream  (directory: IStringStream)



## Test-ITstream  (directory: ITstream)



## Test-IjkField  (directory: IjkField)

    Functionality of IjkField


## Test-IndirectList  (directory: IndirectList)



## Test-IntRange  (directory: IntRange)

    Test integer range

## Test-LabelledItem  (directory: LabelledItem)

    Test LabelledItem (formerly 'Keyed', but that was never used)


## Test-List  (directory: List)

    Simple tests and examples of use of List

See also
    Foam::List


## Test-List2  (directory: List2)

    Test speeds, usability of some List/FixedList operations


## Test-List3  (directory: List3)

    Test list construction


## Test-ListOps  (directory: ListOps)



## Test-ListOps2  (directory: ListOps2)



## Test-ListRead1  (directory: ListRead1)

    List reading


## Test-Map  (directory: Map)



## Test-MathFunctions  (directory: MathFunctions)

    Tests for \c Math namespace member functions
    using \c doubleScalar base type.


## Test-NamedEnum  (directory: NamedEnum)

    Testing of NamedEnum.
    The class is deprecated, but we may still need to support it for things
    like swak4Foam etc.


## Test-OCharStream1  (directory: OCharStream1)



## Test-OCountStream  (directory: OCountStream)

    Test null and counting output streams


## Test-ODE  (directory: ODE)



## Test-OFstream  (directory: OFstream)

    Test OFstream. Primarily atomic operations


## Test-OSspecific  (directory: OSspecific)

    Report some basic os-specific values


## Test-OStringStream  (directory: OStringStream)



## Test-OTstream  (directory: OTstream)



## Test-PDRblockMesh  (directory: PDRblockMesh)

    Test accessors for PDRblock


## Test-PackedList  (directory: PackedList)



## Test-PackedList1  (directory: PackedList1)



## Test-PackedList2  (directory: PackedList2)



## Test-PatchEdgeFaceWave  (directory: PatchEdgeFaceWave)

    Test PatchEdgeFaceWave.


## Test-PatchFunction1  (directory: PatchFunction1)

    Tests Function1


## Test-PatchTools  (directory: PatchTools)

    Test app for PatchTools functionality


## Test-PointEdgeWave  (directory: PointEdgeWave)

    Test pointEdgeWave.


## Test-Polynomial  (directory: Polynomial)

    Test application for the templated Polynomial class


## Test-PrecisionAdaptor  (directory: PrecisionAdaptor)



## Test-PtrList  (directory: PtrList)

    Test behaviour of UPtrList, PtrList

## Test-PtrListDictionary  (directory: PtrListDictionary)



## Test-PtrMap  (directory: PtrMap)



## Test-Random  (directory: Random)

    Simple test for sequence of random numbers


## Test-SLList  (directory: SLList)



## Test-SpanStream1  (directory: SpanStream1)



## Test-SphericalTensor  (directory: SphericalTensor)

    Tests for \c SphericalTensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-SphericalTensor2D  (directory: SphericalTensor2D)

    Tests for \c SphericalTensor2D constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-SubField  (directory: SubField)

    Simple tests on SubList, SubField


## Test-SymmTensor  (directory: SymmTensor)

    Tests for \c SymmTensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c symmTensor, i.e. SymmTensor<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-SymmTensor2D  (directory: SymmTensor2D)

    Tests for \c SymmTensor2D constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c symmTensor2D, i.e. SymmTensor2D<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-Tensor  (directory: Tensor)

    Tests for \c Tensor constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c tensor, i.e. Tensor<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-Tensor2D  (directory: Tensor2D)

    Tests for \c Tensor2D constructors, member functions and operators
    using \c floatScalar, \c doubleScalar, and \c complex base types.

    Eigen decomposition tests for \c tensor2D, i.e. Tensor2D<scalar>.

    Cross-checks were obtained from 'NumPy 1.15.1' and 'SciPy 1.1.0' if no
    theoretical cross-check exists (like eigendecomposition relations), and
    were hard-coded for elementwise comparisons.

    For \c complex base type, the cross-checks do only involve zero imag part.


## Test-Tuple2  (directory: Tuple2)

    Test construction, comparison etc for Tuple2 and Pair.


## Test-UDictionary  (directory: UDictionary)



## Test-UIndirectList  (directory: UIndirectList)



## Test-UList  (directory: UList)

    Simple tests for UList constructors

See also
    Foam::List


## Test-UniformField  (directory: UniformField)

    Test uniform list/field

## Test-argList  (directory: argList)



## Test-autoPtr  (directory: autoPtr)

- no description

## Test-barycentric  (directory: barycentric)

    Some simple tests for barycentric coordinates and transforms


## Test-base64Encoding  (directory: base64)

    Test base64 encoding layer.

    Target values generated with "base64 encode ..." in Google.
    A simple independent source for comparison.


## Test-bitSet1  (directory: bitSet1)

    Basic bitSet characteristics


## Test-bitSet2  (directory: bitSet2)

    Test bitSet functionality


## Test-bitops  (directory: bitops)

    Test some bit-operations.


## Test-boolList  (directory: boolList)

    Test specialized boolList functionality


## Test-boolVector  (directory: boolVector)

    Some simple tests for boolVector


## Test-boundBox  (directory: boundBox)

    Test bounding box behaviour


## Test-boundBox2  (directory: boundBox2)

    Test bounding box behaviour


## Test-broadcastCopy  (directory: broadcastCopy)

    Test file/directory broadcasting


## Test-callback  (directory: callback)



## Test-cellModels  (directory: cellModels)

    Print information about known cellModels


## Test-charList  (directory: charList)

    Some test of UList, List for characters


## Test-checkDecomposePar  (directory: checkDecomposePar)

    Check decomposition from kaffpa (KaHIP) output.
    foamToMetisGraph was likely used for producing the kaffpa input.


## Test-clock  (directory: clock)

    Test some clock-related routines


## Test-codeStream  (directory: codeStream)



## Test-colourTables  (directory: colourTables)

- no description

## Test-complex  (directory: complex)

    Tests for complex numbers


## Test-compoundToken1  (directory: compoundToken1)

    Test token construct assign etc.

## Test-constantFields  (directory: constantFields)

    Simple compilation tests for constant fields


## Test-contiguous  (directory: contiguous)

    Simple test of contiguous data


## Test-coordinateSystem  (directory: coordinateSystem)

    Expand coordinate system definitions


## Test-copyFile  (directory: copyFile)

    Test atomic copyFile as per timeActivatedFileUpdate


## Test-cpluplus1  (directory: cplusplus1)

    Test miscellaneous C++ templates/functionality.


## Test-cpuInfo  (directory: cpuInfo)



## Test-cstring  (directory: cstring)

    Test some string functionality


## Test-cyclic  (directory: cyclic)

    Incompressible CFD code


## Test-decomposedBlockData  (directory: decomposedBlockData)

    Convert decomposedBlockData into its components.


## Test-delete  (directory: delete)



## Test-dictionary  (directory: dictionary)



## Test-dictionary2  (directory: dictionary2)

    Test dictionary insertion and some reading functionality.


## Test-dictionary3  (directory: dictionary3)

    Test expressions and re-expansions


## Test-dictionary4  (directory: dictionary4)

    Test expansion

## Test-dictionaryCopy  (directory: dictionaryCopy)

    Test copying a dictionary with filtering


## Test-dictionaryTokens  (directory: dictionaryTokens)


    Test dictionaryTokens


## Test-dimField  (directory: dimField)

    Simple tests for DimensionedField


## Test-dimensionSet  (directory: dimensionSet)

    Print values of predefined dimensionSets, and some other tests


## Test-dimensionedType  (directory: dimensionedType)

- no description

## Test-dynamicIndexedOctree  (directory: dynamicIndexedOctree)

    Test the construction, insertion and removal of points from the dynamic
    indexed octree.


## Test-dynamicLibrary  (directory: dynamicLibrary)

    Test loading/unloading of libraries


## Test-edges  (directory: edges)

    Simple tests for edges


## Test-ensightFile  (directory: ensightFile)

    check cleanup of ensight file and variable names


## Test-error  (directory: error)



## Test-etcFiles  (directory: etcFiles)

    Test etcFiles functionality.
    Similar to foamEtcFile script, but automatically prunes nonexistent
    directories from the list.


## Test-exprEntry  (directory: exprEntry)

    Read in the given dictionaries and attempt to use exprEntry expansion
    on any strings.

Note
   Since this is only for testing purposes, only handles simple dictionary
   entries without attempting to descend into sub-dicts.


## Test-exprTraits  (directory: exprTraits)

    Basic tests of expression traits


## Test-exprValue  (directory: exprValue)

    Test low-level polymorphic value container (exprValue)


## Test-ExtendedStencil  (directory: extendedStencil)

    Test app for determining extended stencil.


## Test-ExtendedStencil2  (directory: extendedStencil)

    Test app for determining extended cell-to-cell stencil.


## Test-externalFileCoupler  (directory: externalFileCoupler)

    Test of master/slave communication etc.


## Test-faceHashing  (directory: faceHashing)

    Basic tests of face/triFace hashing


## Test-faces  (directory: faces)

    Simple tests for various faces


## Test-fft  (directory: fft)

    Very simple fft tests


## Test-field1  (directory: field1)

    Simple field tests

    Test use of Kahan/Neumaier to extend precision for when running SPDP
    mode. Conclusion is that it is easier/quicker to run these summation
    loops as double precision (i.e. solveScalar).


## Test-fieldDependency  (directory: fieldDependency)

    Test field dependencies.


## Test-fieldMapping  (directory: fieldMapping)

    Test app for mapping of fields.


## Test-fieldTypes  (directory: fieldTypes)

    Print fieldTypes


## Test-fileHandler-dummy  (directory: fileHandler-dummy)

    Simple test of dummy fileOperation


## Test-fileHandler-ranks1  (directory: fileHandler-ranks1)

    Test IO ranks and ranks selection


## Test-fileHandler-writing  (directory: fileHandler-writing)

    Simple test of file writing, including timings


## Test-fileName  (directory: fileName)

    Test some basic fileName functionality


## Test-fileNameClean  (directory: fileNameClean)



## Test-fileNameOS  (directory: fileNameOS)

    Test fileName behaviour, potential OS capabilities etc.

    In the distant future could possibly replace parts with C++ filesystem


## Test-fileOperation1  (directory: fileOperation1)

   Test string parsing and other bits for fileOperation


## Test-findCell-octree  (directory: findCell-octree)

- no description

## Test-findSphereFeatureEdges-octree  (directory: findSphereFeatureEdges-octree)

- no description

## Test-findTimes  (directory: findTimes)



## Test-flatOuput1  (directory: flatOutput1)

    Simple test of FlatOutput


## Test-foamEnv  (directory: foamEnv)

    Test etcFiles functionality.
    Similar to foamEtcFile script, but automatically prunes nonexistent
    directories from the list.


## Test-foamVersion  (directory: foamVersion)

    Print the OpenFOAM version information.


## Test-fstreamPointer  (directory: fstreamPointer)

    Low-level fstream tests


## Test-fvSolutionCombine  (directory: fvSolutionCombine)

    Simple utility for combining fvSolution solution entries.


## Test-fvc  (directory: fvc)

    Finite volume method test code.


## Test-fvc2D  (directory: fvc2D)

    Finite volume method test code for 2-D space.


## Test-gatherValues1  (directory: gatherValues1)

    Test list gather functionality


## Test-globalIndex  (directory: globalIndex)

    Simple tests for the globalIndex class.


## Test-globalMeshData  (directory: globalMeshData)

    Test global point communication


## Test-graph  (directory: graph)

    Test program for making graphs


## Test-graphXi  (directory: graphXi)

    Test program for making graphs


## Test-gravityMeshObject  (directory: gravityMeshObject)

    Test loading of different gravity items


## Test-hashedWordList  (directory: hashedWordList)



## Test-hexRef8  (directory: hexRef8)

    Test app for refinement and unrefinement. Runs a few iterations refining
    and unrefining.


## Test-instant  (directory: instant)

    Test instant, fileNameInstant


## Test-invTensor  (directory: invTensor)

    Tests for regular and corner cases of tensor inversion.


## Test-io  (directory: io)

    Test basic stream functionality


## Test-labelRanges  (directory: labelRanges)

    Test label ranges

## Test-leastSquareGrad  (directory: leastSquareGrad)



## Test-limits  (directory: limits)

    Print some numerical limits.


## Test-liquid  (directory: liquid)



## Test-mapDistributePolyMesh  (directory: mapDistributePolyMesh)

    Test for procAddressing


## Test-MappedPatch  (directory: mappedPatch)

    Test mapped b.c. by mapping face centres (mesh.C().boundaryField()).


## Test-maxMem  (directory: maxMem)

- no description

## Test-memInfo  (directory: memInfo)



## Test-mesh  (directory: mesh)

- no description

## Test-minMax1  (directory: minMax1)

    Test minMax


## Test-minMax2  (directory: minMax2)

    Test-minMax2


## Test-mkdir  (directory: mkdir)



## Test-momentOfInertia  (directory: momentOfInertia)

    Calculates the inertia tensor and principal axes and moments of a
    test face, tetrahedron and cell.


## Test-multiDimPolyFitter  (directory: multiDimPolyFitter)



## Test-mvBak  (directory: mvBak)



## Test-namedDictionary  (directory: namedDictionary)

    Test handling of keyType/dictionary


## Test-nullObject  (directory: nullObject)

    Tests of nullObject


## Test-objectRegistry  (directory: objectRegistry)

    Simple test of objectRegistry functionality.
    Particular focus on the behaviour of subRegistry.


## Test-objectRegistry2  (directory: objectRegistry2)

    Print objectRegistry information, with some additional tests.


## Test-pTraits  (directory: pTraits)



## Test-parallel-broadcast  (directory: parallel-broadcast)

    Test for various broadcast routines.


## Test-parallel-chunks  (directory: parallel-chunks)

    Test for sending contiguous data in chunk-wise.
    Largely mirrors Pstream::exchange or vice versa


## Test-parallel-comm0  (directory: parallel-comm0)

    Very basic checks on standard communicators


## Test-parallel-comm1  (directory: parallel-comm1)

    Checks communication using user-defined communicators


## Test-parallel-comm2  (directory: parallel-comm2)

    Basic communicator tests


## Test-parallel-comm3a  (directory: parallel-comm3a)

    Basic communicator tests


## Test-parallel-comm3b  (directory: parallel-comm3b)

    Basic communicator tests


## Test-parallel-external-init  (directory: parallel-external-init)

    Simulate starting MPI outside of OpenFOAM


## Test-parallel-nbx2  (directory: parallel-nbx2)

    Test for send/receive data


## Test-parallel-nonBlocking  (directory: parallel-nonBlocking)

    Test for various non-blocking parallel routines.


## Test-parallel-waitSome  (directory: parallel-waitSome)

    Test polling versus wait-all for processing receive data.
    Will not see much difference between -wait-all and -no-polling though
    since the master doesn't have enough other work.


## Test-parallel  (directory: parallel)

    Test for various parallel routines.


## Test-passiveParticle  (directory: passiveParticle)

    Test cloud of passive particles.


## Test-patchRegion  (directory: patchRegion)

    Detect point pinches


## Test-plotFunction1  (directory: plotFunction1)

    Plot scalar Function1 entries


## Test-PointField  (directory: pointField)

    For each time calculate the magnitude of velocity.


## Test-polyMeshGeom-speed1  (directory: polyMeshGeom-speed1)

    Simple timing tests for some polyMesh primitives


## Test-predicates  (directory: predicates)

    Simple tests using predicates


## Test-prefixOSstream  (directory: prefixOSstream)



## Test-PrimitivePatch  (directory: primitivePatch)

    Test new primitive patches.


## Test-primitives  (directory: primitives)

    Parsing etc for primitives.


## Test-processorTopology  (directory: processorTopology)

    Test/output processor topology


## Test-quaternion  (directory: quaternion)

    Test application for quaternions.


## Test-rawIOField  (directory: rawIOField)

    Reading rawIOField from disk


## Test-readBroadcast1  (directory: readBroadcast1)

    Test file reading with broadcast


## Test-readDir  (directory: readDir)

    Test functionality of Foam::readDir


## Test-reconstruct  (directory: reconstruct)

- no description

## Test-reconstructedDistanceFunction  (directory: reconstructedDistanceFunction)



## Test-refPtr  (directory: refPtr)

    Tests some basic functionality of refPtr


## Test-regex1  (directory: regex1)

    Tests for regular expressions


## Test-processorRouter  (directory: router)



## Test-scalarOps  (directory: scalarOps)

    Test scalar-only ops


## Test-scalarPredicates  (directory: scalarPredicates)

    Simple tests using predicates for scalars


## Test-scalarRanges  (directory: scalarRanges)

    Test scalar ranges

## Test-searchableSphere  (directory: searchableSphere)

    Basic tests for searchable sphere


## Test-SHA1  (directory: sha1)



## Test-sigFpe  (directory: sigFpe)

    Test handling of floating point exceptions by provoking them


## Test-simpleMatrix  (directory: simpleMatrix)

- no description

## Test-sizeof  (directory: sizeof)

    Test the sizeof various classes.


## Test-sliceRange  (directory: sliceRange)

    Test slice range

## Test-slicedField  (directory: slicedField)



## Test-sortList  (directory: sort)



## Test-stdFoam-span  (directory: span)

    Basic functionality test for span


## Test-spline  (directory: spline)

- no description

## Test-splitFunctionArgs  (directory: splitFunctionArgs)

    Test splitting of function name args


## Test-string  (directory: string)

    Test some string functionality


## Test-string2  (directory: string2)

    Test some string functionality


## Test-stringList  (directory: stringList)



## Test-stringSplit  (directory: stringSplit)

    Test string splitting


## Test-string_view1  (directory: string_view1)

    Test some string_view functionality


## Test-surface-sampling  (directory: surface-sampling)

    Simple test of surface sampling, including timings


## Test-surfaceIntersection  (directory: surfaceIntersection)

    Test surface-surface intersection


## Test-surfaceMeshConvert  (directory: surfaceMeshConvert)

    Test conversions from one surface mesh format to another.

Usage
    \b Test-surfaceMeshConvert inputFile outputFile [OPTION]

    Options:
      - \par -clean
        Perform some surface checking/cleanup on the input surface

      - \par -orient
        Check face orientation on the input surface

      - \par -testModify
        Test modification mechanism

      - \par -scale \<scale\>
        Specify a scaling factor for writing the files

      - \par -triSurface
        Use triSurface library for input/output
Note
    The filename extensions are used to determine the file format type.


## Test-surfaceReading  (directory: surfaceReading)

    Test basic surface format reading capabilities (and speeds)

Note
    The filename extensions are used to determine the file format type.


## Test-surfaceTree  (directory: surfaceTree)

    Simple tests for building indexedOctree etc.


## Test-surfaceWriter  (directory: surfaceWriter)

    Test surface writers.

Usage
    \b Test-surfaceWriter inputFile outputFile


## Test-syncTools  (directory: syncTools)

    Test some functionality in syncTools.


## Test-sysInfo  (directory: sysInfo)



## Test-tensor2D  (directory: tensor2D)

    Tests for tensor2D and vector2D


## Test-tensorFields1  (directory: tensorFields1)

- no description

## Test-tetTetOverlap  (directory: tetTetOverlap)

    Overlap volume of two tets


## Test-thermoMixture  (directory: thermoMixture)



## Test-timeSelector  (directory: timeSelector)

    Test TimePaths and timeSelector


## Test-tmp  (directory: tmp)

    Tests for possible memory leaks in the tmp (and tmp<Field> algebra).


## Test-token  (directory: token)

    Test token construct assign etc.

## Test-tokenize  (directory: tokenize)

    Test the tokenizing of various things

## Test-treeComms  (directory: treeComms)

    Print/test tree communication patterns


## Test-triTet  (directory: triTet)


## Test-triangleIntersection  (directory: triangleIntersection)

    Test bounding box / triangle intersection


## Test-unitConversion  (directory: unitConversion)

- no description

## Test-vector  (directory: vector)

    Some simple tests for vector


## Test-vectorTools  (directory: vectorTools)


## Test-volField  (directory: volField)

- no description

## Test-volPointInterpolation  (directory: volPointInterpolation)

- no description

## Test-vtkSeriesWriter  (directory: vtkSeriesWriter)

    Basic functionality tests for vtk::seriesWriter


## Test-vtmWriter  (directory: vtmWriter)

    Basic functionality tests for vtk::vtmWriter


## Test-wallDist  (directory: wallDist)

    Calculate and write the distance-to-wall field for a moving mesh.


## Test-wallDistDyM  (directory: wallDistDyM)

    Calculate and write the distance-to-wall field for a moving mesh.


## Test-wmake1  (directory: wmake1)

    Some tests for wmake features.
    For example, testing how robust or fragile version-dependent conditional
    compilation works.


## Test-wordRe  (directory: wordRe)

    Test word/regex


## Test-write-wrapped-string  (directory: write-wrapped-string)

    Simple tests for wrapped strings


## Test-zoneDistribute  (directory: zoneDistribute)

    Test of zoneDistribute validated with mapDistribute

    Original code supplied by Henning Scheufler, DLR (2019)


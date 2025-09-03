/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define with the build version of dill */
#define DILL_VERSION "3.0.0"

/* Define if bfd functions use bfd_byte arguments. */
#undef BFD_BYTE

/* Define to 1 if the emulator should be built */
/* #undef BUILD_EMULATOR */

/* Define if dill should attempt to DCG */
/* #undef EMULATION_ONLY */

/* Define to 1 if you have the <dis-asm.h> header file. */
/* #undef HAVE_DIS_ASM_H */

/* Define to 1 if you have the <dlfcn.h> header file. */
/* #undef HAVE_DLFCN_H */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* Define to 1 if you have __clear_cache is defined  */
#define CLEAR_CACHE_DEFINED

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H

/* Define to 1 if you have the <sys/mman.h> header file. */
#define HAVE_SYS_MMAN_H

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H

/* Define if you have the `print_insn_arm' function. */
/* #undef HAVE_PRINT_INSN_ARM */

/* Define if you have the `print_insn_big_powerpc' function. */
/* #undef HAVE_PRINT_INSN_BIG_POWERPC */

/* Define if you have the `print_insn_i386' function. */
/* #undef HAVE_PRINT_INSN_I386 */

/* Define if you have the `print_insn_ia64' function. */
/* #undef HAVE_PRINT_INSN_IA64 */

/* Define if you have the `print_insn_little_arm' function. */
/* #undef HAVE_PRINT_INSN_LITTLE_ARM */

/* Define if you have the `print_insn_sparc' function. */
/* #undef HAVE_PRINT_INSN_SPARC */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H

/* Define to 1 if you have the <sys/stat.h> header file. */
/* #undef HAVE_SYS_STAT_H */

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H

/* Define if the host processor is ARM5 */
/* #undef HOST_ARM5 */

/* Define if the host processor is ARM6 */
/* #undef HOST_ARM6 */

/* Define if the host processor is ARM7 */
/* #undef HOST_ARM7 */

/* Define if the host processor is ARM8 */
/* #undef HOST_ARM8 */

/* Define if the host processor is an ia64 */
/* #undef HOST_IA64 */

/* Define if the host processor is a powerpc */
/* #undef HOST_POWERPC */

/* Define if the host processor is a powerpc64 */
/* #undef HOST_POWERPC64 */

/* Define if the host processor is a ppc64le */
/* #undef HOST_PPC64LE */

/* Define if the host processor is a sparc */
/* #undef HOST_SPARC */

/* Define if the host processor is a sparcv9 */
/* #undef HOST_SPARCV9 */

/* Define if the host processor is an x86 */
/* #undef HOST_X86 */

/* Define if the host processor is an x86_64 */
#define HOST_X86_64

/* Define if INIT_DISASSEMBLE_INFO takes three arguments instead of two */
/* #undef INIT_DISASSEMBLE_INFO_THREE_ARG */

/* Define if integrating with kernel plugins */
/* #undef KPLUGINS_INTEGRATION */

/* Define if compiling for linux kernel */
/* #undef LINUX_KERNEL_MODULE */

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
/* #undef LT_OBJDIR */

/* Define if you want more than just native target */
/* #undef MULTI_TARGET */

/* Define for the host architecture type */
#define NATIVE_ARCH "x86_64"

/* Define if there is no disassembler */
#define NO_DISASSEMBLER

/* Define if we should not use inlined procedures from BFD */
/* #undef NO_INLINED_BFD_PROCS */

/* The number of bytes in type long */
#define SIZEOF_LONG 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS

/* Define if the membar instruction should be used to sync Icache and Dcache
   */
/* #undef USE_MEMBAR */

/* Define this if mmap should be used instead of malloc() for code memory */
#define USE_MMAP_CODE_SEG

/* Define this if VirtualProtect should be used to change memory protections */
/* #undef USE_VIRTUAL_PROTECT */

/* Define this if windows calling convention should be used */
/* #undef USE_WINDOWS_CALLS */

/* Define if byteorder is bigendian */
/* #undef WORDS_BIGENDIAN */

/* Define if using ARM hardware float ABIs for parameter passing */
#define ARM_HARD_FLOAT  0

/* Define if using ARMv7 hardware */
/* #undef ARMV7_AVAILABLE */

/* Define if we shouldn't use native DCG */
/* #undef DILL_IGNORE_NATIVE */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
#undef inline
#endif

// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause

/*

This file provides a unified front-end for the wrapper generators.

*/

#include "vtkParseMain.h"
#include "vtkParse.h"
#include "vtkParseData.h"
#include "vtkParseDependencyTracking.h"
#include "vtkParseSystem.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

/* This is the struct that contains the options */
OptionInfo options;

/* Get the base filename */
static const char* parse_exename(const char* cmd)
{
  const char* exename;

  /* remove directory part of exe name */
  for (exename = cmd + strlen(cmd); exename > cmd; --exename)
  {
    char pc = exename[-1];
    if (pc == ':' || pc == '/' || pc == '\\')
    {
      break;
    }
  }

  return exename;
}

/* Print the help */
static void parse_print_help(FILE* fp, const char* cmd, int multi)
{
  fprintf(fp,
    "Usage: %s [options] infile... \n"
    "  --help            print this help message\n"
    "  --version         print the VTK version\n"
    "  -o <file>         the output file\n"
    "  -I <dir>          add an include directory\n"
    "  -D <macro[=def]>  define a preprocessor macro\n"
    "  -U <macro>        undefine a preprocessor macro\n"
    "  -imacros <file>   read macros from a header file\n"
    "  -MF <file>        write dependency information to a file\n"
    "  -undef            do not predefine platform macros\n"
    "  -Wempty           warn when nothing is wrapped\n"
    "  -Wno-empty        do not warn when nothing is wrapped\n"
    "  @<file>           read arguments from a file\n",
    parse_exename(cmd));

  /* args for describing a single header file input */
  if (!multi)
  {
    fprintf(fp,
      "  -dM               dump all macro definitions to output\n"
      "  --hints <file>    the hints file to use\n"
      "  --types <file>    the type hierarchy file to use\n");
  }
}

/* append an arg to the arglist */
static void parse_append_arg(int* argn, char*** args, char* arg)
{
  /* if argn is a power of two, allocate more space */
  if (*argn > 0 && (*argn & (*argn - 1)) == 0)
  {
    *args = (char**)realloc(*args, 2 * (*argn) * sizeof(char*));
  }
  /* append argument to list */
  (*args)[*argn] = arg;
  (*argn)++;
}

/* read options from a file, return zero on error */
static int read_option_file(StringCache* strings, const char* filename, int* argn, char*** args)
{
  static int option_file_stack_max = 10;
  static int option_file_stack_size = 0;
  static const char* option_file_stack[10];
  FILE* fp;
  char* line;
  const char* ccp;
  char* argstring;
  char* arg;
  size_t maxlen = 15;
  size_t i, n;
  int j;
  int in_string;

  /* TODO: track this dependency properly; tracking is never active at this
   * point. */
  fp = vtkParse_FileOpen(filename, "r");

  if (fp == NULL)
  {
    return 0;
  }

  line = (char*)malloc(maxlen);

  /* read the file line by line */
  while (fgets(line, (int)maxlen, fp))
  {
    n = strlen(line);

    /* if buffer not long enough, increase it */
    while (n == maxlen - 1 && line[n - 1] != '\n' && !feof(fp))
    {
      char* oldline = line;
      maxlen *= 2;
      line = (char*)realloc(line, maxlen);
      if (!line)
      {
        free(oldline);
        fclose(fp);
        return 0;
      }
      if (!fgets(&line[n], (int)(maxlen - n), fp))
      {
        break;
      }
      n += strlen(&line[n]);
    }

    /* allocate a string to hold the parsed arguments */
    argstring = vtkParse_NewString(strings, n);
    arg = argstring;
    i = 0;

    /* break the line into individual options */
    ccp = line;
    in_string = 0;
    while (*ccp != '\0')
    {
      for (;;)
      {
        if (*ccp == '\\')
        {
          ccp++;
        }
        else if (*ccp == '\"' || *ccp == '\'')
        {
          if (!in_string)
          {
            in_string = *ccp++;
            continue;
          }
          else if (*ccp == in_string)
          {
            in_string = 0;
            ccp++;
            continue;
          }
        }
        else if (!in_string && isspace(*ccp))
        {
          do
          {
            ccp++;
          } while (isspace(*ccp));
          break;
        }
        if (*ccp == '\0')
        {
          break;
        }
        /* append character to argument */
        arg[i++] = *ccp++;
      }
      arg[i++] = '\0';

      if (arg[0] == '@')
      {
        /* recursively expand '@file' option */
        if (option_file_stack_size == option_file_stack_max)
        {
          fprintf(stderr, "%s: @file recursion is too deep.\n", (*args)[0]);
          exit(1);
        }
        /* avoid reading the same file recursively */
        option_file_stack[option_file_stack_size++] = filename;
        for (j = 0; j < option_file_stack_size; j++)
        {
          if (strcmp(&arg[1], option_file_stack[j]) == 0)
          {
            break;
          }
        }
        if (j < option_file_stack_size)
        {
          parse_append_arg(argn, args, arg);
        }
        else if (read_option_file(strings, &arg[1], argn, args) == 0)
        {
          parse_append_arg(argn, args, arg);
        }
        option_file_stack_size--;
      }
      else if (arg[0] != '\0')
      {
        parse_append_arg(argn, args, arg);
      }
      /* prepare for next arg */
      arg += i;
      i = 0;
    }
  }

  free(line);
  fclose(fp);
  return 1;
}

/* expand any "@file" args that occur in the command-line args */
static void parse_expand_args(StringCache* strings, int argc, char* argv[], int* argn, char*** args)
{
  int i;

  *argn = 0;
  *args = (char**)malloc(sizeof(char*));

  for (i = 0; i < argc; i++)
  {
    /* check for "@file" unless this is the command name */
    if (i > 0 && argv[i][0] == '@')
    {
      /* if read_option_file returns null, add "@file" to the args */
      /* (this mimics the way that gcc expands @file arguments) */
      if (read_option_file(strings, &argv[i][1], argn, args) == 0)
      {
        parse_append_arg(argn, args, argv[i]);
      }
    }
    else
    {
      /* append any other arg */
      parse_append_arg(argn, args, argv[i]);
    }
  }
}

/* Check the options: "multi" should be zero for wrapper tools that
 * only take one input file, or one for wrapper tools that take multiple
 * input files.  Returns zero for "--version" or "--help", or returns -1
 * if an error occurred.  Otherwise, it returns the number of args
 * that were successfully parsed. */
static int parse_check_options(int argc, char* argv[], int multi)
{
  int i;
  size_t j;
  char* cp;
  char c;

  options.NumberOfFiles = 0;
  options.Files = NULL;
  options.InputFileName = NULL;
  options.OutputFileName = NULL;
  options.NumberOfHierarchyFileNames = 0;
  options.HierarchyFileNames = NULL;
  options.NumberOfHintFileNames = 0;
  options.HintFileNames = NULL;
  options.DumpMacros = 0;
  options.DepFileName = NULL;
  options.WarningFlags.Empty = 0;

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i], "--help") == 0)
    {
      parse_print_help(stdout, argv[0], multi);
      return 0;
    }
    else if (strcmp(argv[i], "--version") == 0)
    {
      const char* ver = VTK_PARSE_VERSION;
      fprintf(stdout, "%s %s\n", parse_exename(argv[0]), ver);
      return 0;
    }
    else if (strcmp(argv[i], "-imacros") == 0)
    {
      i++;
      if (i >= argc || argv[i][0] == '-')
      {
        return -1;
      }
      cp = argv[i];
      vtkParse_IncludeMacros(cp);
    }
    else if (strcmp(argv[i], "-undef") == 0)
    {
      vtkParse_UndefinePlatformMacros();
    }
    else if (strcmp(argv[i], "-dM") == 0)
    {
      options.DumpMacros = 1;
    }
    else if (strcmp(argv[i], "-Wempty") == 0)
    {
      options.WarningFlags.Empty = 1;
    }
    else if (strcmp(argv[i], "-Wno-empty") == 0)
    {
      options.WarningFlags.Empty = 0;
    }
    else if (strcmp(argv[i], "-MF") == 0)
    {
      i++;
      if (i >= argc || argv[i][0] == '-')
      {
        return -1;
      }
      options.DepFileName = argv[i];
    }
    else if (argv[i][0] == '-' && isalpha(argv[i][1]))
    {
      c = argv[i][1];
      cp = &argv[i][2];
      if (*cp == '\0')
      {
        i++;
        if (i >= argc || argv[i][0] == '-')
        {
          return -1;
        }
        cp = argv[i];
      }

      if (c == 'o')
      {
        options.OutputFileName = cp;
      }
      else if (c == 'I')
      {
        vtkParse_IncludeDirectory(cp);
      }
      else if (c == 'D')
      {
        j = 0;
        while (cp[j] != '\0' && cp[j] != '=')
        {
          j++;
        }
        if (cp[j] == '=')
        {
          j++;
        }
        vtkParse_DefineMacro(cp, &cp[j]);
      }
      else if (c == 'U')
      {
        vtkParse_UndefineMacro(cp);
      }
    }
    else if (argv[i][0] != '-')
    {
      if (options.NumberOfFiles == 0)
      {
        options.Files = (char**)malloc(sizeof(char*));
      }
      else if ((options.NumberOfFiles & (options.NumberOfFiles - 1)) == 0)
      {
        options.Files = (char**)realloc(options.Files, 2 * options.NumberOfFiles * sizeof(char*));
      }
      options.Files[options.NumberOfFiles++] = argv[i];
    }
    else if (!multi && strcmp(argv[i], "--hints") == 0)
    {
      i++;
      if (i >= argc || argv[i][0] == '-')
      {
        return -1;
      }
      if (options.NumberOfHintFileNames == 0)
      {
        options.HintFileNames = (char**)malloc(sizeof(char*));
      }
      else if ((options.NumberOfHintFileNames & (options.NumberOfHintFileNames - 1)) == 0)
      {
        options.HintFileNames =
          (char**)realloc(options.HintFileNames, 2 * options.NumberOfHintFileNames * sizeof(char*));
      }
      options.HintFileNames[options.NumberOfHintFileNames++] = argv[i];
    }
    else if (!multi && strcmp(argv[i], "--types") == 0)
    {
      i++;
      if (i >= argc || argv[i][0] == '-')
      {
        return -1;
      }
      if (options.NumberOfHierarchyFileNames == 0)
      {
        options.HierarchyFileNames = (char**)malloc(sizeof(char*));
      }
      else if ((options.NumberOfHierarchyFileNames & (options.NumberOfHierarchyFileNames - 1)) == 0)
      {
        options.HierarchyFileNames = (char**)realloc(
          options.HierarchyFileNames, 2 * options.NumberOfHierarchyFileNames * sizeof(char*));
      }
      options.HierarchyFileNames[options.NumberOfHierarchyFileNames++] = argv[i];
    }
  }

  return i;
}

/* Return a pointer to the static OptionInfo struct */
const OptionInfo* vtkParse_GetCommandLineOptions(void)
{
  return &options;
}

int vtkParse_Finalize(void)
{
  int ret = 0;

  if (options.DepFileName && vtkParse_DependencyTrackingWrite(options.DepFileName))
  {
    ret = 1;
  }
  vtkParse_FinalizeDependencyTracking();

  return ret;
}

/* Command-line argument handler for wrapper tools */
FileInfo* vtkParse_Main(int argc, char* argv[])
{
  int argi;
  FILE* ifile;
  FILE* hfile = 0;
  int nhfiles;
  int ihfiles;
  const char* hfilename;
  FileInfo* data;
  StringCache strings;
  int argn;
  char** args;

  /* set the command name for diagnostics */
  vtkParse_SetCommandName(parse_exename(argv[0]));

  /* hook the cleanup function */
  atexit(vtkParse_FinalCleanup);

  /* pre-define the __VTK_WRAP__ macro */
  vtkParse_DefineMacro("__VTK_WRAP__", 0);

  /* expand any "@file" args */
  vtkParse_InitStringCache(&strings);
  parse_expand_args(&strings, argc, argv, &argn, &args);

  /* read the args into the static OptionInfo struct */
  argi = parse_check_options(argn, args, 0);

  /* verify number of args, print usage if not valid */
  if (argi == 0)
  {
    free(args);
    exit(0);
  }
  else if (argi < 0 || options.NumberOfFiles != 1)
  {
    parse_print_help(stderr, args[0], 0);
    exit(1);
  }

  /* open the input file */
  options.InputFileName = options.Files[0];

  if (!(ifile = vtkParse_FileOpen(options.InputFileName, "r")))
  {
    fprintf(stderr, "Error opening input file %s\n", options.InputFileName);
    exit(1);
  }

  /* free the expanded args */
  free(args);

  /* make sure that an output file was given on the command line,
   * unless dumping info in which case stdout will be used instead */
  if (options.DumpMacros)
  {
    vtkParse_DumpMacros(options.OutputFileName);
  }
  else if (options.OutputFileName == NULL)
  {
    fprintf(stderr, "No output file was specified\n");
    fclose(ifile);
    exit(1);
  }

  if (options.DepFileName && options.OutputFileName)
  {
    vtkParse_InitDependencyTracking(options.OutputFileName);
    /* TODO: register response files read in `read_option_file` here. */
  }

  /* parse the input file */
  data = vtkParse_ParseFile(options.InputFileName, ifile, stderr);

  if (!data)
  {
    vtkParse_FreeStringCache(&strings);
    exit(1);
  }

  /* merge into a single string cache to avoid leaking strings */
  vtkParse_MergeStringCache(data->Strings, &strings);

  /* check whether -dM option was set */
  if (options.DumpMacros)
  {
    /* do nothing (the dump occurred in ParseFile above) */
    exit(0);
  }

  /* open and parse each hint file, if given on the command line */
  nhfiles = options.NumberOfHintFileNames;
  for (ihfiles = 0; ihfiles < nhfiles; ihfiles++)
  {
    hfilename = options.HintFileNames[ihfiles];
    if (hfilename && hfilename[0] != '\0')
    {
      if (!(hfile = vtkParse_FileOpen(hfilename, "r")))
      {
        fprintf(stderr, "Error opening hint file %s\n", hfilename);
        fclose(ifile);
        vtkParse_FreeFile(data);
        exit(1);
      }

      /* fill in some blanks by using the hints file */
      vtkParse_ReadHints(data, hfile, stderr);
    }
  }

  if (data->MainClass)
  {
    /* mark class as abstract unless it has New() method */
    int nfunc = data->MainClass->NumberOfFunctions;
    int ifunc;
    for (ifunc = 0; ifunc < nfunc; ifunc++)
    {
      const FunctionInfo* func = data->MainClass->Functions[ifunc];
      if (func && func->Access == VTK_ACCESS_PUBLIC && func->Name &&
        strcmp(func->Name, "New") == 0 && func->NumberOfParameters == 0)
      {
        break;
      }
    }
    data->MainClass->IsAbstract = ((ifunc == nfunc) ? 1 : 0);
  }

  return data;
}

/* Command-line argument handler for wrapper tools */
StringCache* vtkParse_MainMulti(int argc, char* argv[])
{
  int argi;
  int argn;
  char** args;
  StringCache* strings = (StringCache*)malloc(sizeof(StringCache));

  /* set the command name for diagnostics */
  vtkParse_SetCommandName(parse_exename(argv[0]));

  /* hook the cleanup function */
  atexit(vtkParse_FinalCleanup);

  /* pre-define the __VTK_WRAP__ macro */
  vtkParse_DefineMacro("__VTK_WRAP__", 0);

  /* expand any "@file" args */
  vtkParse_InitStringCache(strings);
  parse_expand_args(strings, argc, argv, &argn, &args);

  /* read the args into the static OptionInfo struct */
  argi = parse_check_options(argn, args, 1);
  free(args);

  if (argi == 0)
  {
    exit(0);
  }
  else if (argi < 0 || options.NumberOfFiles == 0)
  {
    parse_print_help(stderr, argv[0], 1);
    exit(1);
  }

  if (options.DepFileName && options.OutputFileName)
  {
    vtkParse_InitDependencyTracking(options.OutputFileName);
    /* TODO: register response files read in `read_option_file` here. */
  }

  /* the input file */
  options.InputFileName = options.Files[0];
  return strings;
}

#ifdef _WIN32

/* To hold the wmain() args after conversion to narrow args */
static char** parse_win32_argv = NULL;

/* Cleanup function to be called at exit, after wmain() */
static void parse_win32_cleanup(void)
{
  free(parse_win32_argv);
}

/* Convert wchar_t args to utf8 */
char** vtkParse_WideArgsToUTF8(int argc, wchar_t* wargv[])
{
  int i, n;
  int cl = 0;
  char** argv;
  char* cp;

  /* compute total command-line length */
  for (i = 0; i < argc; i++)
  {
    cl += WideCharToMultiByte(CP_UTF8, 0, wargv[i], -1, NULL, 0, NULL, NULL);
  }

  /* allocate combined buffer for argv and arg strings */
  argv = (char**)malloc(argc * sizeof(char*) + cl);
  cp = (char*)(argv + argc);

  /* convert all arguments */
  for (i = 0; i < argc; i++)
  {
    argv[i] = cp;
    n = WideCharToMultiByte(CP_UTF8, 0, wargv[i], -1, argv[i], cl, NULL, NULL);
    cp += n;
    cl -= n;
  }

  /* use atexit() to handle the cleanup */
  parse_win32_argv = argv;
  atexit(parse_win32_cleanup);

  return argv;
}
#endif

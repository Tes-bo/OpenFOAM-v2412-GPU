set(_VisItReader_cmake_dir "${CMAKE_CURRENT_LIST_DIR}")

function (add_visit_reader NAME VERSION)
  set(_add_visit_reader_plugin_name "vtk${NAME}")

  cmake_parse_arguments(_add_visit_reader
    ""
    "VISIT_READER_NAME;VISIT_INCLUDE_NAME;VISIT_READER_TYPE;VISIT_READER_USES_OPTIONS;VISIT_READER_OPTIONS_NAME"
    "SERVER_SOURCES"
    ${ARGN})

  if (_add_visit_reader_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR
      "Unparsed arguments for add_visit_reader: "
      "${_add_visit_reader_UNPARSED_ARGUMENTS}")
  endif ()

  if (NOT DEFINED _add_visit_reader_VISIT_READER_NAME)
    message(FATAL_ERROR
      "The `VISIT_READER_NAME` argument is required.")
  endif ()

  if (NOT DEFINED _add_visit_reader_VISIT_READER_USES_OPTIONS)
    set(_add_visit_reader_VISIT_READER_USES_OPTIONS OFF)
  endif ()

  # If the user hasn't defined an include name, we presume the reader name is
  # also the include name.
  if (NOT DEFINED _add_visit_reader_VISIT_INCLUDE_NAME)
    set(_add_visit_reader_VISIT_INCLUDE_NAME "${_add_visit_reader_VISIT_READER_NAME}")
  endif ()

  if (NOT _add_visit_reader_VISIT_READER_TYPE MATCHES "^[SM]T[SM]D$")
    message(FATAL_ERROR
      "Invalid Reader Type `${_add_visit_reader_VISIT_READER_TYPE}`. Valid "
      "Types are STSD, STMD, MTSD, MTMD")
  endif ()

  if (_add_visit_reader_VISIT_READER_USES_OPTIONS)
    # Determine the name of the plugin info class by removing the AVT from the
    # start and the FileFormat from the back
    if (_add_visit_reader_VISIT_READER_OPTIONS_NAME)
      set(_add_visit_reader_name "${_add_visit_reader_VISIT_READER_OPTIONS_NAME}")
    else ()
      set(_add_visit_reader_name "${_add_visit_reader_VISIT_READER_NAME}")
    endif ()
    string(REGEX REPLACE "^avt|FileFormat$" "" _add_visit_reader_name_base "${_add_visit_reader_name}")
    set(_add_visit_reader_visit_plugin_info_header
      "${_add_visit_reader_name_base}PluginInfo")
    set(_add_visit_reader_visit_plugin_info_class
      "${_add_visit_reader_name_base}CommonPluginInfo")
  endif ()

  # Need to generate the VTK class wrapper.
  string(SUBSTRING "${_add_visit_reader_VISIT_READER_TYPE}" 0 2 _add_visit_reader_wrapper_type)

  set(_add_visit_reader_source
    "${_VisItReader_output_dir}/${_add_visit_reader_plugin_name}.cxx")
  set(_add_visit_reader_header
    "${_VisItReader_output_dir}/${_add_visit_reader_plugin_name}.h")

  configure_file(
    "${_VisItReader_cmake_dir}/VisIt${_add_visit_reader_wrapper_type}.h.in"
    "${_add_visit_reader_header}"
    @ONLY)
  configure_file(
    "${_VisItReader_cmake_dir}/VisIt${_add_visit_reader_wrapper_type}.cxx.in"
    "${_add_visit_reader_source}"
    @ONLY)

  # Fix up the arg_server_sources path for compilation
  set(_add_visit_reader_server_sources)
  foreach (_add_visit_reader_server_source IN LISTS _add_visit_reader_SERVER_SOURCES)
    if (NOT IS_ABSOLUTE "${_add_visit_reader_server_source}")
      string(PREPEND _add_visit_reader_server_source
        "${CMAKE_CURRENT_SOURCE_DIR}/")
    endif ()
    list(APPEND _add_visit_reader_server_sources
      "${_add_visit_reader_server_source}")
  endforeach ()

  set_property(GLOBAL APPEND
    PROPERTY
      _visit_server_sources "${_add_visit_reader_server_sources}")
  set_property(GLOBAL APPEND
    PROPERTY
      _visit_db_sources "${_add_visit_reader_source}")
  set_property(GLOBAL APPEND
    PROPERTY
      _visit_db_headers "${_add_visit_reader_header}")
  set_property(GLOBAL APPEND
    PROPERTY
      _visit_db_include_dirs "${CMAKE_CURRENT_SOURCE_DIR}")
endfunction ()

function(ADD_VISIT_INTERFACE_READER NAME VERSION)
  set(INTERFACE_NAME "vtk${NAME}")
  set(INTERFACE_VERSION "${VERSION}")

  cmake_parse_arguments(_add_visit_reader
    ""
    "VISIT_INTERFACE_CALL;VISIT_INTERFACE_FILE"
    "VISIT_READERS;SERVER_SOURCES"
    ""
    ${ARGN})

  if (_add_visit_reader_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR
      "Unparsed arguments for add_visit_interface_reader: "
      "${_add_visit_reader_UNPARSED_ARGUMENTS}")
  endif ()

  if (NOT DEFINED _add_visit_reader_VISIT_INTERFACE_CALL)
    message(FATAL_ERROR
      "The `VISIT_INTERFACE_CALL` argument is required.")
  endif ()

  if (NOT DEFINED _add_visit_reader_VISIT_INTERFACE_FILE)
    message(FATAL_ERROR
      "The `VISIT_INTERFACE_FILE` argument is required.")
  endif ()

  set(_add_visit_reader_parse "")
  set(_add_visit_reader_names)
  set(_add_visit_reader_name)
  set(_add_visit_reader_type)
  set(_add_visit_reader_include)
  foreach (_add_visit_reader_arg IN LISTS _add_visit_reader_VISIT_READERS)
    if (_add_visit_reader_parse STREQUAL "")
      set(_add_visit_reader_parse "${_add_visit_reader_arg}")
    elseif (_add_visit_reader_parse STREQUAL "NAME")
      set(_add_visit_reader_name "${_add_visit_reader_arg}")
      list(APPEND _add_visit_reader_names "${_add_visit_reader_name}")
      set(_add_visit_reader_parse "")
    elseif (_add_visit_reader_parse STREQUAL "USES_INTERFACE")
      if (NOT _add_visit_reader_name)
        message(FATAL_ERROR
          "Missing `NAME`")
      endif ()
      if (DEFINED "_add_visit_reader_${_add_visit_reader_name}_uses_interface")
        message(FATAL_ERROR
          "Duplicate `USES_INTERFACE` for `${_add_visit_reader_name}`")
      endif ()
      set("_add_visit_reader_${_add_visit_reader_name}_uses_interface"
        "${_add_visit_reader_arg}")
      set(_add_visit_reader_parse "")
    elseif (_add_visit_reader_parse STREQUAL "TYPE")
      if (NOT _add_visit_reader_name)
        message(FATAL_ERROR
          "Missing `NAME`")
      endif ()
      if (DEFINED "_add_visit_reader_${_add_visit_reader_name}_type")
        message(FATAL_ERROR
          "Duplicate `TYPE` for `${_add_visit_reader_name}`")
      endif ()
      set("_add_visit_reader_${_add_visit_reader_name}_type"
        "${_add_visit_reader_arg}")
      set(_add_visit_reader_parse "")
    elseif (_add_visit_reader_parse STREQUAL "INCLUDE")
      if (NOT _add_visit_reader_name)
        message(FATAL_ERROR
          "Missing `NAME`")
      endif ()
      if (DEFINED "_add_visit_reader_${_add_visit_reader_name}_include")
        message(FATAL_ERROR
          "Duplicate `INCLUDE` for `${_add_visit_reader_name}`")
      endif ()
      set("_add_visit_reader_${_add_visit_reader_name}_include"
        "${_add_visit_reader_arg}")
      set(_add_visit_reader_parse "")
    else ()
      message(FATAL_ERROR
        "Unknown argument `${_add_visit_reader_parse}`")
    endif ()
  endforeach ()

  if (_add_visit_reader_parse)
    message(FATAL_ERROR
      "Missing argument for `${_add_visit_reader_parse}`")
  endif ()

  if (NOT _add_visit_reader_names)
    message(FATAL_ERROR
      "No `NAME` arguments given")
  endif ()

  set(_add_visit_reader_sources)
  set(_add_visit_reader_headers)
  foreach (_add_visit_reader_name IN LISTS _add_visit_reader_names)
    if (NOT DEFINED "_add_visit_reader_${_add_visit_reader_name}_type")
      message(FATAL_ERROR
        "Missing `TYPE` for `${_add_visit_reader_name}`")
    endif ()
    if (NOT DEFINED "_add_visit_reader_${_add_visit_reader_name}_uses_interface")
      set("_add_visit_reader_${_add_visit_reader_name}_uses_interface" ON)
    endif ()
    if (NOT DEFINED "_add_visit_reader_${_add_visit_reader_name}_include")
      set("_add_visit_reader_${_add_visit_reader_name}_include"
        "${_add_visit_reader_name}")
    endif ()

    set(_add_visit_reader_VISIT_READER_NAME
      "${_add_visit_reader_name}")
    set(_add_visit_reader_type
      "${_add_visit_reader_${_add_visit_reader_name}_type}")
    set(_add_visit_reader_uses_interface
      "${_add_visit_reader_${_add_visit_reader_name}_uses_interface}")
    set(_add_visit_reader_VISIT_INCLUDE_NAME
      "${_add_visit_reader_${_add_visit_reader_name}_include}")

    string(REGEX REPLACE "^avt|FileFormat$" "" _add_visit_reader_name_base "${_add_visit_reader_name}")
    set(_add_visit_reader_plugin_name "vtkVisIt${_add_visit_reader_name_base}Reader")

    string(SUBSTRING "${_add_visit_reader_type}" 0 2 _add_visit_reader_wrapper_type)

    set(_add_visit_reader_interface_header
      "${_VisItReader_output_dir}/${_add_visit_reader_plugin_name}${_add_visit_reader_VISIT_INTERFACE_FILE}.h")
    set(_add_visit_reader_source
      "${_VisItReader_output_dir}/${_add_visit_reader_plugin_name}.cxx")
    set(_add_visit_reader_header
      "${_VisItReader_output_dir}/${_add_visit_reader_plugin_name}.h")

    configure_file(
      "${CMAKE_CURRENT_SOURCE_DIR}/${_add_visit_reader_VISIT_INTERFACE_FILE}.h.in"
      "${_add_visit_reader_interface_header}"
      @ONLY)
    configure_file(
      "${_VisItReader_cmake_dir}/VisIt${_add_visit_reader_wrapper_type}.h.in"
      "${_add_visit_reader_header}"
      @ONLY)
    configure_file(
      "${_VisItReader_cmake_dir}/VisIt${_add_visit_reader_wrapper_type}.cxx.in"
      "${_add_visit_reader_source}"
      @ONLY)

    list(APPEND _add_visit_reader_sources
      "${_add_visit_reader_interface_header}"
      "${_add_visit_reader_source}")
    list(APPEND _add_visit_reader_headers
      "${_add_visit_reader_header}")
  endforeach ()

  # Fix up the arg_server_sources path for compilation
  set(_add_visit_reader_server_sources)
  foreach (_add_visit_reader_server_source IN LISTS _add_visit_reader_SERVER_SOURCES)
    if (NOT IS_ABSOLUTE "${_add_visit_reader_server_source}")
      string(PREPEND _add_visit_reader_server_source
        "${CMAKE_CURRENT_SOURCE_DIR}/")
    endif ()
    list(APPEND _add_visit_reader_server_sources
      "${_add_visit_reader_server_source}")
  endforeach ()

  set_property(GLOBAL APPEND
    PROPERTY
      _visit_server_sources "${_add_visit_reader_server_sources}")
  set_property(GLOBAL APPEND
    PROPERTY
      _visit_db_sources "${_add_visit_reader_sources}")
  set_property(GLOBAL APPEND
    PROPERTY
      _visit_db_headers "${_add_visit_reader_headers}")
  set_property(GLOBAL APPEND
    PROPERTY
      _visit_db_include_dirs "${CMAKE_CURRENT_SOURCE_DIR}")
endfunction ()

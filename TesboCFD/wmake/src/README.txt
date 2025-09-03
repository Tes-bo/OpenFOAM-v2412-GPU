Minimal changes to lemon for C++ integration.

- Additional '-e' command line option to define the code extension.
  By default this is 'c', but with this option can define something
  like '-ecxx' etc for using a C++ compiler.

- Additional '%static' Lemon directive, which is boolean-like:

      %static

  This adds a 'static' qualifier to all of the 'Parse*' routines that
  would otherwise have global linkage, thus making them only visible
  in the same file-scope.
  Can subsequently place all of the C++ interface code within a %code
  block in the same translation unit.

  This allows good encapsulation without fundamentally changing how
  Lemon works.

--
2020-07-10

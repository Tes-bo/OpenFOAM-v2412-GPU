# VTK-m data Directory #

  + This directory contains data and baseline images for VTK-m regression testing.
  + The data/data directory are data files of various types. This includes
    rectilinear and unstructured data.
  + data/data/sentinel-data is a file that is used by CMake to verify if a
    target machine has git lfs installed. If git lfs is not installed,
    then testing will be automatically disabled.
  + The data/baseline are the testing images. These are used in testing to
    compare a valid image against a generated image. If a difference between
    the two images is found, then the test is considered to have failed.

## adiosDict specification

* To simplify the support of multi-regions, multiple clouds, multiple
  species, the default entries of the `adiosDict` can be used.

---

#### The `write` sub-dictionary block

The adios output treats an empty or missing `write` sub-dictionary
block as equivalent to the following example and should work directly
without any user intervention:

    // Additional items to write, or items to write (when explicit)
    write
    {
        // Implicitly selected fields/clouds (where AUTO_WRITE is on)
        explicit    false;

        // Expand cloud components (debugging)
        expandCloud false;

        // No additional volume fields
        fields      ();

        // No additional clouds
        clouds      ();
    }


#### The `ignore` sub-dictionary block

The `ignore` sub-dictionary block can be used to blacklist automatically
(and even explicitly) selected fields/clouds. Blacklisting entire
regions is currently not supported, but could be a possible future
extension.

    // Items to suppress from writing
    ignore
    {
        fields  ();
        clouds  ();
        // Future? -> regions ();
    }


#### The `regions` sub-dictionary block

The `regions` sub-dictionary block can be used to provide fine-grained
control of individual regions. It can currently contain both `write`
and `ignore` sub-blocks. If these are not present, the values of the
enclosing dictionary are used.

This could be used, for example, to block writing of particular
fields in a region.

    regions
    {
        solidRegion
        {
            ignore  { fields (".*"); }
        }
    }

Note that this mechanism cannot currently be used to suppress writing
the mesh for any particular region.

#### Adios write methods

    // Output write method: 'adios_config -m' to see available methods
    // = MPI               (default)
    // = MPI_AGGREGATE
    // = POSIX
    writeMethod     MPI;

    // Parameters for the selected write method
    // = ""               (default)
    // = "verbose=4;num_aggregators=2;num_ost=2" ...
    writeOptions    "";


#### Adios read methods

Not fully implemented

    // Input read method:
    //  = BP               (default)
    //  = BP_AGGREGATE
    readMethod      BP;


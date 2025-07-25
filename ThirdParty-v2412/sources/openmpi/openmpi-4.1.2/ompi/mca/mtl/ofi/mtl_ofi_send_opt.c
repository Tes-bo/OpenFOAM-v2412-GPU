/*
 * Copyright (c) 2013-2018 Intel, Inc. All rights reserved
 *
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "mtl_ofi.h"

__opal_attribute_always_inline__ static inline int
ompi_mtl_ofi_send_false(struct mca_mtl_base_module_t *mtl,
                  struct ompi_communicator_t *comm,
                  int dest,
                  int tag,
                  struct opal_convertor_t *convertor,
                  mca_pml_base_send_mode_t mode)
{
    const bool OFI_CQ_DATA = false;

    return ompi_mtl_ofi_send_generic(mtl, comm, dest, tag,
                                    convertor, mode,
                                    OFI_CQ_DATA);
}

__opal_attribute_always_inline__ static inline int
ompi_mtl_ofi_send_true(struct mca_mtl_base_module_t *mtl,
                  struct ompi_communicator_t *comm,
                  int dest,
                  int tag,
                  struct opal_convertor_t *convertor,
                  mca_pml_base_send_mode_t mode)
{
    const bool OFI_CQ_DATA = true;

    return ompi_mtl_ofi_send_generic(mtl, comm, dest, tag,
                                    convertor, mode,
                                    OFI_CQ_DATA);
}

void ompi_mtl_ofi_send_symtable_init(struct ompi_mtl_ofi_symtable *sym_table)
{

    sym_table->ompi_mtl_ofi_send[false]
        = ompi_mtl_ofi_send_false;


    sym_table->ompi_mtl_ofi_send[true]
        = ompi_mtl_ofi_send_true;

}


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
ompi_mtl_ofi_iprobe_false(struct mca_mtl_base_module_t *mtl,
                struct ompi_communicator_t *comm,
                int src,
                int tag,
                int *flag,
                struct ompi_status_public_t *status)
{
    const bool OFI_CQ_DATA = false;

    return ompi_mtl_ofi_iprobe_generic(mtl, comm, src, tag,
                                    flag, status,
                                    OFI_CQ_DATA);
}

__opal_attribute_always_inline__ static inline int
ompi_mtl_ofi_iprobe_true(struct mca_mtl_base_module_t *mtl,
                struct ompi_communicator_t *comm,
                int src,
                int tag,
                int *flag,
                struct ompi_status_public_t *status)
{
    const bool OFI_CQ_DATA = true;

    return ompi_mtl_ofi_iprobe_generic(mtl, comm, src, tag,
                                    flag, status,
                                    OFI_CQ_DATA);
}

void ompi_mtl_ofi_iprobe_symtable_init(struct ompi_mtl_ofi_symtable *sym_table)
{

    sym_table->ompi_mtl_ofi_iprobe[false]
        = ompi_mtl_ofi_iprobe_false;


    sym_table->ompi_mtl_ofi_iprobe[true]
        = ompi_mtl_ofi_iprobe_true;

}


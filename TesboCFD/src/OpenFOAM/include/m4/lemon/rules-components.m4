divert(-1)dnl
#-----------------------------------*- m4 -*-----------------------------------
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | www.openfoam.com
#      \\/     M anipulation  |
#------------------------------------------------------------------------------
#     Copyright (C) 2019 OpenCFD Ltd.
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM, distributed under GPL-3.0-or-later.
#
# Description
#     Collection of VectorSpace `component' functions
#
#     `rule_scalar_components'
#     `rule_vector_components'
#     `rule_sphTensor_components'
#     `rule_symTensor_components'
#     `rule_tensor_components'
#     `rule_tensor_unzipAll'  (diag and rows)
#
# Defined after inclusion
#     `rule_method_component'
#     `rule_tensor_unzipDiag'
#     `rule_tensor_unzipRow'
#
#------------------------------------------------------------------------------

# These are to be defined *after* inclusion - undefine now

undefine([rule_method_component])
undefine([rule_tensor_unzipDiag])
undefine([rule_tensor_unzipRow])


#------------------------------------------------------------------------------
# rule_scalar_components(source)
#
# Description
#     Extract scalar field from scalar field - no-op
#------------------------------------------------------------------------------

define([rule_scalar_components], [])


#------------------------------------------------------------------------------
# rule_vector_components(out, in)
#
# Description
#     Extract scalar field from vector field
#------------------------------------------------------------------------------

define([rule_vector_components],
[rule_method_component($1, $2, CMPT_X, Foam::vector::X)
rule_method_component($1, $2, CMPT_Y, Foam::vector::Y)
rule_method_component($1, $2, CMPT_Z, Foam::vector::Z)]
)


#------------------------------------------------------------------------------
# rule_sphTensor_components(out, in)
#
# Description
#     Extract scalar field from sphericalTensor field
#------------------------------------------------------------------------------

define([rule_sphTensor_components],
[rule_method_component($1, $2, CMPT_II, Foam::sphericalTensor::II)]
)


#------------------------------------------------------------------------------
# rule_symTensor_components(out, in)
#
# Description
#     Extract scalar field from symmTensor field
#------------------------------------------------------------------------------

define([rule_symTensor_components],
[rule_method_component($1, $2, CMPT_XX, Foam::symmTensor::XX)
rule_method_component($1, $2, CMPT_XY, Foam::symmTensor::XY)
rule_method_component($1, $2, CMPT_XZ, Foam::symmTensor::XZ)
rule_method_component($1, $2, CMPT_YY, Foam::symmTensor::YY)
rule_method_component($1, $2, CMPT_YZ, Foam::symmTensor::YZ)
rule_method_component($1, $2, CMPT_ZZ, Foam::symmTensor::ZZ)]
)

#------------------------------------------------------------------------------
# rule_tensor_components(out, in)
#
# Description
#     Extract scalar field from tensor field
#------------------------------------------------------------------------------

define([rule_tensor_components],
[rule_method_component($1, $2, CMPT_XX, Foam::tensor::XX)
rule_method_component($1, $2, CMPT_XY, Foam::tensor::XY)
rule_method_component($1, $2, CMPT_XZ, Foam::tensor::XZ)
rule_method_component($1, $2, CMPT_YX, Foam::tensor::YX)
rule_method_component($1, $2, CMPT_YY, Foam::tensor::YY)
rule_method_component($1, $2, CMPT_YZ, Foam::tensor::YZ)
rule_method_component($1, $2, CMPT_ZX, Foam::tensor::ZX)
rule_method_component($1, $2, CMPT_ZY, Foam::tensor::ZY)
rule_method_component($1, $2, CMPT_ZZ, Foam::tensor::ZZ)]
)


#------------------------------------------------------------------------------
# rule_tensor_unzipAll(out, in)
#
# Description
#     Extract vector diagonal and rows from tensor field
#------------------------------------------------------------------------------

define([rule_tensor_unzipAll],
[rule_tensor_unzipDiag($1, $2)
rule_tensor_unzipRow($1, $2, CMPT_X, Foam::vector::X)
rule_tensor_unzipRow($1, $2, CMPT_Y, Foam::vector::Y)
rule_tensor_unzipRow($1, $2, CMPT_Z, Foam::vector::Z)]
)


#------------------------------------------------------------------------------
divert(0)dnl

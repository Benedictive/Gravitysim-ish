set(HEADER_LIST "${GravitySim_SOURCE_DIR}/include/gravitysim/kernel.hpp")

add_library(gravitysim_library kernel.cu ${HEADER_LIST})

target_include_directories(gravitysim_library PUBLIC ../include)

set_target_properties(gravitysim_library PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(gravitysim_library PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
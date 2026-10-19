set(SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF
    torsion_constants.c
    protocol_constants.c
)

add_library(${LIB_PRECOMP_${SVARIANT_UPPER}} ${SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE common ${PROJECT_SOURCE_DIR}/protocols/ref/include ${PROJECT_SOURCE_DIR}/protocols/soa/include ${PROJECT_SOURCE_DIR}/ec/ref/include ${PROJECT_SOURCE_DIR}/ec/soa/include ${PROJECT_SOURCE_DIR}/ec/ref/${SVARIANT_LOWER}/include ${PROJECT_SOURCE_DIR}/gf/ref/${SVARIANT_LOWER}/include ${PROJECT_SOURCE_DIR}/gf/soa/${SVARIANT_LOWER}/include ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_COMMON})
target_compile_options(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
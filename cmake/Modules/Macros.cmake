
function(GET_KMER_DEFS VAR)
    math(EXPR MAX_KMER_SIZE "((${HIPMER_KMER_LEN} + 31) / 32) * 32")
    math(EXPR MAX_KMER_PACKED_LENGTH "(${MAX_KMER_SIZE} + 3) / 4")
    set(_TMP MAX_KMER_SIZE=${MAX_KMER_SIZE} MAX_KMER_PACKED_LENGTH=${MAX_KMER_PACKED_LENGTH})
    set(${VAR} ${_TMP} PARENT_SCOPE)
endfunction()

#moved here from main cmake by me
macro(ADD_READID_DEFS_TO_LIB TARGET)
    set(READID_DEFS -DMAX_NUM_READS=${MAX_NUM_READS})
    target_compile_definitions(${TARGET} PUBLIC ${READID_DEFS})
endmacro()

macro(ADD_READID_DEFS TARGET)
    add_readid_defs_to_lib(${TARGET})
endmacro()
# int8 	Byte (-128 to 127)
# int16 	Integer (-32768 to 32767)
# int32 	Integer (-2147483648 to 2147483647)
# int64 	Integer (-9223372036854775808 to 9223372036854775807)
# uint8 	Unsigned integer (0 to 255)
# uint16 	Unsigned integer (0 to 65535)
# uint32 	Unsigned integer (0 to 4294967295)
# uint64 	Unsigned integer (0 to 18446744073709551615)
GENOMIC_COORD_DTYPE = "uint32"  # should be fine as long as individual chromosomes are less than 4Gb
READ_COORD_DTYPE = "uint32"
STRAND_DTYPE = "bool"
FRAG_IDX_DTYPE = "uint32"
READ_IDX_DTYPE = "uint32"
ALIGN_IDX_DTYPE = "uint32"
PERCENTAGE_DTYPE = "float32"
HAPLOTYPE_IDX_DTYPE = "int8"
MQ_DTYPE = "uint8"
ALIGN_SCORE_DTYPE = "uint32"
PHASE_SET_DTYPE = GENOMIC_COORD_DTYPE
HAPLOTYPE_IDX_DTYPE = "int8"
PERCENTAGE_DTYPE = "float32"


INPUT_REFGENOME_REGEX = r"(.+)\.(fasta|fa|fna)(\.gz)*"

PQ_ENGINE = "pyarrow"
PQ_VERSION = "2.0"

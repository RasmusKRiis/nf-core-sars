/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BUILD_PRIMER_DB } from '../modules/local/build_primer_db/main'
include { PRIMER_MISMATCH } from '../modules/local/primer_mismatch/main'
include { SPLIT_FASTA     } from '../modules/local/split_fasta/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parsePrimerFastaInputs(fastaPattern) {
    Channel
        .fromPath(fastaPattern, checkIfExists: true)
        .ifEmpty { error "No FASTA files found with pattern: ${fastaPattern}" }
        .map { fasta ->
            def sampleId = fasta.getBaseName().replaceAll(/\.(fa|fasta|fna)$/, '')
            tuple([id: sampleId], fasta)
        }
}

def resolveRequiredFile(pathValue, label) {
    if (!pathValue) {
        error "${label} was not provided."
    }
    def resolved = file(pathValue)
    if (!resolved.exists()) {
        error "${label} not found: ${pathValue}"
    }
    return resolved
}

def toPathList(paths) {
    paths instanceof List ? paths : [paths]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW DEFINITION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SARSCOVSEQPRIMERCHECK {

    main:
    def fastaPattern = params.fasta ?: error('Please provide --fasta pointing to the input multi-FASTA file(s).')
    def primerBedFile = resolveRequiredFile(params.primer_bed, 'Primer BED file')
    def primerFastaFile = resolveRequiredFile(params.primer_fasta, 'Primer FASTA file')
    def primerSetName = params.primer_set_name ?: primerBedFile.getBaseName()

    ch_fasta_raw = parsePrimerFastaInputs(fastaPattern)

    SPLIT_FASTA(ch_fasta_raw)

    ch_single_seq = SPLIT_FASTA.out.split_fastas
        .flatMap { meta, splitFastas ->
            toPathList(splitFastas).collect { splitFasta ->
                tuple([id: splitFasta.getBaseName()], splitFasta)
            }
        }

    BUILD_PRIMER_DB(
        Channel.value(primerBedFile),
        Channel.value(primerFastaFile),
        Channel.value(primerSetName)
    )

    PRIMER_MISMATCH(
        ch_single_seq,
        BUILD_PRIMER_DB.out.primer_db,
        Channel.value(params.runid)
    )
}

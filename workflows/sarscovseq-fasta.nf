/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CSV_CONVERSION } from '../modules/local/csv_conversion/main'
include { NEXTCLADE      } from '../modules/local/nextclade/main'
include { REPORTFASTA    } from '../modules/local/report-fasta/main'
include { SPLIT_FASTA    } from '../modules/local/split_fasta/main'
include { TABLELOOKUP    } from '../modules/local/tablelookup/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def parseFastaInputs(fastaPattern) {
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

workflow SARSCOVSEQFASTA {

    main:
    def fastaPattern = params.fasta ?: error('Please provide --fasta pointing to one or more FASTA files.')
    def spikeTable = resolveRequiredFile(params.spike, 'Spike resistance lookup table')
    def rdrpTable = resolveRequiredFile(params.rdrp, 'RdRp resistance lookup table')
    def clproTable = resolveRequiredFile(params.clpro, '3CLpro resistance lookup table')

    ch_fasta_raw = parseFastaInputs(fastaPattern)

    SPLIT_FASTA(ch_fasta_raw)

    ch_single_seq = SPLIT_FASTA.out.split_fastas
        .flatMap { meta, splitFastas ->
            toPathList(splitFastas).collect { splitFasta ->
                tuple([id: splitFasta.getBaseName()], splitFasta)
            }
        }

    NEXTCLADE(ch_single_seq)

    CSV_CONVERSION(NEXTCLADE.out.nextclade_csv)

    TABLELOOKUP(
        CSV_CONVERSION.out.nextclade_mutations,
        Channel.value(spikeTable),
        Channel.value(rdrpTable),
        Channel.value(clproTable)
    )

    REPORTFASTA(
        CSV_CONVERSION.out.nextclade_stats_report.collect(),
        CSV_CONVERSION.out.nextclade_mutations_report.collect(),
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        params.runid
    )
}

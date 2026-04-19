/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_FASTQ          } from '../modules/nf-core/cat/fastq/main'
include { ARTIC_GUPPYPLEX    } from '../modules/local/artic_guppyplex/main'
include { ARTIC_MINION_M     } from '../modules/local/artic_minion_m/main'
include { BUILD_PRIMER_DB    } from '../modules/local/build_primer_db/main'
include { CHOPPER            } from '../modules/local/chopper/main'
include { CSV_CONVERSION     } from '../modules/local/csv_conversion/main'
include { DEPTH_ANALYSIS     } from '../modules/local/depth_analysis/main'
include { NEXTCLADE          } from '../modules/local/nextclade/main'
include { PRIMER_MISMATCH    } from '../modules/local/primer_mismatch/main'
include { REPORT             } from '../modules/local/report/main'
include { TABLELOOKUP        } from '../modules/local/tablelookup/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

def collectSampleFastqs(samplesDirPath, barcode) {
    def barcodeDir = new File(samplesDirPath.toString(), barcode)
    if (!barcodeDir.exists() || !barcodeDir.isDirectory()) {
        return []
    }

    return (barcodeDir.listFiles() ?: [])
        .findAll { candidate -> candidate.name.endsWith('.fastq.gz') || candidate.name.endsWith('.fastq') }
        .sort { left, right -> left.name <=> right.name }
        .collect { candidate -> file(candidate.toString()) }
}

def parseSampleSheet(sampleSheetPath) {
    def sampleSheetFile = file(sampleSheetPath)
    if (!sampleSheetFile.exists()) {
        error "Input samplesheet not found: ${sampleSheetPath}"
    }

    def samplesDirPath = file(params.samplesDir)
    if (!samplesDirPath.exists()) {
        error "The samples directory defined by --samplesDir does not exist: ${samplesDirPath}"
    }
    if (!samplesDirPath.toFile().isDirectory()) {
        error "The path defined by --samplesDir is not a directory: ${samplesDirPath}"
    }

    return Channel
        .fromPath(sampleSheetFile.toString(), checkIfExists: true)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            if (!row.SequenceID || !row.Barcode) {
                error "Missing 'SequenceID' or 'Barcode' in the samplesheet row: ${row}"
            }

            def sampleId = row.SequenceID.toString().trim()
            def barcode = row.Barcode.toString().trim()
            def fastqs = collectSampleFastqs(samplesDirPath, barcode)

            if (!fastqs) {
                log.warn "Skipping sample '${sampleId}': no FASTQ files found in ${samplesDirPath}/${barcode}"
                return null
            }

            def meta = [id: sampleId, barcode: barcode, single_end: true]
            tuple(meta, fastqs)
        }
        .filter { it != null }
        .ifEmpty {
            error "No samples from ${sampleSheetFile} matched FASTQ files under ${samplesDirPath}"
        }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SARSCOVSEQ {

    main:
    def sampleSheetFile = resolveRequiredFile(params.input, 'Input samplesheet')
    def spikeTable = resolveRequiredFile(params.spike, 'Spike resistance lookup table')
    def rdrpTable = resolveRequiredFile(params.rdrp, 'RdRp resistance lookup table')
    def clproTable = resolveRequiredFile(params.clpro, '3CLpro resistance lookup table')

    def primerDirPath = params.primerdir ? file(params.primerdir) : null

    def primerBedCandidates = []
    if (params.primer_bed) {
        primerBedCandidates << file(params.primer_bed)
    }
    if (primerDirPath?.exists()) {
        def primerDirFiles = primerDirPath.toFile().listFiles() ?: []
        primerBedCandidates << file("${params.primerdir}/SARS-CoV-2.scheme.bed")
        primerBedCandidates << file("${params.primerdir}/ncov-2019_midnight.scheme.bed")
        primerBedCandidates.addAll(primerDirFiles.findAll { it.name ==~ /.*\.scheme\.bed$/ })
        primerBedCandidates.addAll(primerDirFiles.findAll { it.name ==~ /.*\.bed$/ })
    }
    def primerBedCandidate = primerBedCandidates.find { it && it.exists() }
    def primerBedFile = primerBedCandidate ? file(primerBedCandidate.toString()) : null
    if (!primerBedFile) {
        error "Unable to locate primer BED file. Provide --primer_bed or place a .scheme.bed file inside --primerdir."
    }

    def referenceCandidates = []
    if (params.reference) {
        referenceCandidates << file(params.reference)
    }
    if (primerDirPath?.exists()) {
        def primerDirFiles = primerDirPath.toFile().listFiles() ?: []
        referenceCandidates << file("${params.primerdir}/SARS-CoV-2.reference.fasta")
        referenceCandidates << file("${params.primerdir}/ncov-2019_midnight.reference.fasta")
        referenceCandidates.addAll(primerDirFiles.findAll { it.name ==~ /.*\.reference\.fasta$/ })
    }
    def referenceCandidate = referenceCandidates.find { it && it.exists() }
    def referenceFile = referenceCandidate ? file(referenceCandidate.toString()) : null
    if (!referenceFile) {
        error "Unable to locate reference FASTA. Provide --reference or place a *.reference.fasta inside --primerdir."
    }

    def primerFastaCandidates = []
    if (params.primer_fasta) {
        def explicitPrimerFasta = file(params.primer_fasta)
        if (!explicitPrimerFasta.exists()) {
            error "Primer FASTA file not found: ${params.primer_fasta}"
        }
        primerFastaCandidates << explicitPrimerFasta
    }
    if (primerDirPath?.exists()) {
        primerFastaCandidates << file("${params.primerdir}/primers.fasta")
        primerFastaCandidates << file("${params.primerdir}/SARS-CoV-2.primers.fasta")
    }
    def bedParent = primerBedFile.getParent()
    if (bedParent) {
        primerFastaCandidates << file("${bedParent}/primers.fasta")
        primerFastaCandidates << file("${bedParent}/SARS-CoV-2.primers.fasta")
    }
    def primerFastaFile = primerFastaCandidates.find { it && it.exists() }
    if (!primerFastaFile) {
        error "Unable to locate a primers FASTA file. Provide --primer_fasta or place primers.fasta alongside the selected primer scheme."
    }

    def inferredPrimerName = params.primer_set_name
    if (!inferredPrimerName) {
        if (primerDirPath?.exists()) {
            inferredPrimerName = primerDirPath.getFileName()?.toString()
        }
        if (!inferredPrimerName && bedParent) {
            inferredPrimerName = bedParent.getFileName()?.toString()
        }
        if (!inferredPrimerName) {
            def fileName = primerBedFile.getFileName()?.toString() ?: 'primer_scheme'
            inferredPrimerName = fileName.replaceAll(/\.bed$/, '')
        }
    }
    def primerSetName = (inferredPrimerName ?: 'primer_set').replaceAll(/[^A-Za-z0-9._-]+/, '_')
    def primerDisplayPath = primerDirPath ? primerDirPath.toString() : primerBedFile.toString()

    ch_sample_information = parseSampleSheet(sampleSheetFile)

    BUILD_PRIMER_DB(
        Channel.value(primerBedFile),
        Channel.value(primerFastaFile),
        Channel.value(primerSetName)
    )

    CAT_FASTQ(ch_sample_information)

    CHOPPER(CAT_FASTQ.out.reads)

    ARTIC_GUPPYPLEX(CHOPPER.out.chopperfastq)

    ARTIC_MINION_M(
        ARTIC_GUPPYPLEX.out.gp_fastq,
        Channel.value(primerBedFile),
        Channel.value(referenceFile)
    )

    DEPTH_ANALYSIS(
        ARTIC_MINION_M.out.artic_bam,
        Channel.value(primerBedFile)
    )

    PRIMER_MISMATCH(
        ARTIC_MINION_M.out.artic_consensus,
        BUILD_PRIMER_DB.out.primer_db,
        Channel.value(params.runid)
    )

    NEXTCLADE(ARTIC_MINION_M.out.artic_consensus)

    CSV_CONVERSION(NEXTCLADE.out.nextclade_csv)

    TABLELOOKUP(
        CSV_CONVERSION.out.nextclade_mutations,
        Channel.value(spikeTable),
        Channel.value(rdrpTable),
        Channel.value(clproTable)
    )

    REPORT(
        CSV_CONVERSION.out.nextclade_stats_report.collect(),
        CSV_CONVERSION.out.nextclade_mutations_report.collect(),
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        params.runid,
        params.release_version,
        params.seq_instrument,
        Channel.value(sampleSheetFile),
        primerDisplayPath,
        ARTIC_MINION_M.out.artic_consensus_report.collect()
    )
}

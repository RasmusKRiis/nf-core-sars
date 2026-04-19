
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_sarscovseq_pipeline'

include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { AMPLIGONE                   } from '../modules/local/ampligone/main'
include { CHOPPER                     } from '../modules/local/chopper/main'
include { IRMA                        } from '../modules/local/irma/main'
include { NEXTCLADE                   } from '../modules/local/nextclade/main'
include { CSV_CONVERSION              } from '../modules/local/csv_conversion/main'
include { TABLELOOKUP                 } from '../modules/local/tablelookup/main'
include { REPORT                      } from '../modules/local/report/main'
include { BUILD_PRIMER_DB             } from '../modules/local/build_primer_db/main'
include { PRIMER_MISMATCH             } from '../modules/local/primer_mismatch/main'
include { DEPTH_ANALYSIS              } from '../modules/local/depth_analysis/main'

include { PRIMER_MASK                 } from '../modules/local/primer_mask/main'
include { MK_MASK_NO_PRIMER           } from '../modules/local/mk_mask_no_primer/main'
include { MEDAKA_VARIANT              } from '../modules/local/medaka_variant/main'
include { BCFTOOLS_CONSENSUS          } from '../modules/local/bcftools_consensus/main'

include { ARTIC_GUPPYPLEX             } from '../modules/local/artic_guppyplex/main'
include { ARTIC_MINION_M              } from '../modules/local/artic_minion_m/main'







/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Function to parse the sample sheet
def parseSampleSheet(sampleSheetPath) {
    return Channel
        .fromPath(sampleSheetPath)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            // Use 'SequenceID' and 'Barcode' based on the sample sheet
            if (!row.SequenceID || !row.Barcode) {
                error "Missing 'SequenceID' or 'Barcode' in the sample sheet row: ${row}"
            }

            // Use SequenceID as sample ID and Barcode to find files
            def sampleId = row.SequenceID
            def files = file("${params.samplesDir}/${row.Barcode}/*.fastq.gz")

            // Check if there are any files in the list
            if (!files || files.size() == 0) {
                log.warn "Skipping sample ${sampleId}: No FASTQ files found in ${files}"
                return null // Return null to allow filtering out later
            }

            // Creating a metadata map
            def meta = [ id: sampleId, single_end: true ]
            return tuple(meta, files)
        }
        .filter { it != null } // Remove null entries (samples with no FASTQ files)
}


workflow SARSCOVSEQ() {

    main:

    def primerDirPath = params.primerdir ? file(params.primerdir) : null
    def primerBedCandidates = []
    if( params.primer_bed ) {
        primerBedCandidates << file(params.primer_bed)
    }
    if( primerDirPath?.exists() ) {
        def primerDirFiles = primerDirPath.toFile().listFiles() ?: []
        primerBedCandidates << file("${params.primerdir}/SARS-CoV-2.scheme.bed")
        primerBedCandidates << file("${params.primerdir}/ncov-2019_midnight.scheme.bed")
        primerBedCandidates.addAll(primerDirFiles.findAll { it.name ==~ /.*\\.scheme\\.bed$/ })
        primerBedCandidates.addAll(primerDirFiles.findAll { it.name ==~ /.*\\.bed$/ })
    }
    def primerBedCandidate = primerBedCandidates.find { it && it.exists() }
    def primerBedFile = primerBedCandidate ? file(primerBedCandidate.toString()) : null
    if( !primerBedFile ) {
        error "Unable to locate primer BED file. Provide --primer_bed or place a .scheme.bed file inside --primerdir."
    }
    if( !primerBedFile.exists() ) {
        error "Primer BED file not found: ${primerBedFile}"
    }

    def referenceCandidates = []
    if( params.reference ) {
        referenceCandidates << file(params.reference)
    }
    if( primerDirPath?.exists() ) {
        def primerDirFiles = primerDirPath.toFile().listFiles() ?: []
        referenceCandidates << file("${params.primerdir}/SARS-CoV-2.reference.fasta")
        referenceCandidates << file("${params.primerdir}/ncov-2019_midnight.reference.fasta")
        referenceCandidates.addAll(primerDirFiles.findAll { it.name ==~ /.*\\.reference\\.fasta$/ })
    }
    def referenceCandidate = referenceCandidates.find { it && it.exists() }
    def referenceFile = referenceCandidate ? file(referenceCandidate.toString()) : null
    if( !referenceFile ) {
        error "Unable to locate reference FASTA. Provide --reference or place a *.reference.fasta inside --primerdir."
    }
    if( !referenceFile.exists() ) {
        error "Reference FASTA file not found: ${referenceFile}"
    }

    def primerFastaCandidates = []
    if( params.primer_fasta ) {
        def explicitPrimerFasta = file(params.primer_fasta)
        if( !explicitPrimerFasta.exists() ) {
            error "Primer FASTA file not found: ${params.primer_fasta}"
        }
        primerFastaCandidates << explicitPrimerFasta
    }
    if( primerDirPath?.exists() ) {
        primerFastaCandidates << file("${params.primerdir}/primers.fasta")
        primerFastaCandidates << file("${params.primerdir}/SARS-CoV-2.primers.fasta")
    }
    def bedParent = primerBedFile.getParent()
    if( bedParent ) {
        primerFastaCandidates << file("${bedParent}/primers.fasta")
        primerFastaCandidates << file("${bedParent}/SARS-CoV-2.primers.fasta")
    }
    def primerFastaFile = primerFastaCandidates.find { it && it.exists() }
    if( !primerFastaFile ) {
        error "Unable to locate a primers FASTA file. Provide --primer_fasta or ensure a primers.fasta is placed alongside the selected primer directory."
    }

    def inferredPrimerName = params.primer_set_name
    if( !inferredPrimerName ) {
        if( primerDirPath?.exists() ) {
            inferredPrimerName = primerDirPath.getFileName()?.toString()
        }
        if( !inferredPrimerName && bedParent ) {
            inferredPrimerName = bedParent.getFileName()?.toString()
        }
        if( !inferredPrimerName ) {
            def fileName = primerBedFile.getFileName()?.toString() ?: "primer_scheme"
            inferredPrimerName = fileName.replaceAll(/\\.bed$/, '')
        }
    }
    def primerSetName = (inferredPrimerName ?: "primer_set").replaceAll(/[^A-Za-z0-9._-]+/, "_")

    def currentDir = System.getProperty('user.dir')



    ch_sample_information = parseSampleSheet(params.input) // Use params.input directly
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
 

    ch_sample_information
        .map { meta, files ->
            tuple(meta, files.toList())
        }
    .set { read_input }

    BUILD_PRIMER_DB(
        Channel.value(primerBedFile),
        Channel.value(primerFastaFile),
        Channel.value(primerSetName)
    )


    //
    // MODULE: RUN CAT FASTQ
    //
    CAT_FASTQ (
        read_input
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: CHOPPER
    //  
    
    CHOPPER (
        CAT_FASTQ.out.reads
    )


    // MODULE ARTIC

    ARTIC_GUPPYPLEX(
    CHOPPER.out.chopperfastq
    )

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


    //
    // MODULE: AMPLIGONE. 
    //    

    //AMPLIGONE (
    //    CHOPPER.out.chopperfastq, Channel.value(file(params.primerdir))
    //)


    //
    // MODULE: IRMA
    //

    //IRMA (
    //    AMPLIGONE.out.primertrimmedfastq
    //)


    //
    // MODULE: DEPTH ANALYSIS
    //

    //DEPTH_ANALYSIS (
    //    IRMA.out.bam

    //)


    //
    // MODULE: NEXTCLADE
    //

    NEXTCLADE (
        //IRMA.out.amended_consensus
        ARTIC_MINION_M.out.artic_consensus
    
    )

    //
    // MODULE: NEXTCLADE CONVERSION
    //

    CSV_CONVERSION (
        NEXTCLADE.out.nextclade_csv
    )

    //
    // MODULE: NEXTCLADE CONVERSION
    //

    TABLELOOKUP (
        CSV_CONVERSION.out.nextclade_mutations, Channel.value(file(params.spike)), Channel.value(file(params.rdrp)), Channel.value(file(params.clpro))
    )

    //
    // MODULE: REPORT
    //

    def runid = params.runid
    def seq_instrument   = params.seq_instrument
    def primer   = "${params.primerdir}" 
    def release_version = params.release_version


    REPORT (
        CSV_CONVERSION.out.nextclade_stats_report.collect(), 
        CSV_CONVERSION.out.nextclade_mutations_report.collect(), 
        TABLELOOKUP.out.resistance_mutations_report.collect(),
        runid,
        release_version,
        seq_instrument,
        Channel.value(file(params.input)),
        primer,
        //IRMA.out.amended_consensus_report.collect()
        ARTIC_MINION_M.out.artic_consensus_report.collect()
    )



    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    read_input
    //)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

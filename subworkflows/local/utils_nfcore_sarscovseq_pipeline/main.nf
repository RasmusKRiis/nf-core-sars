//
// Subworkflow with functionality specific to the nf-core/sars pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap } from 'plugin/nf-schema'
include { completionEmail } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version
    validate_params
    monochrome_logs
    nextflow_cli_args
    outdir
    input
    help
    help_full
    show_hidden

    main:
    ch_versions = channel.empty()

    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/sars ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    ${workflow.manifest.homePage}/blob/${workflow.manifest.defaultBranch}/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.homePage} -profile <docker/singularity/.../institute> --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    validateInputParameters()

    emit:
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email
    email_on_fail
    plaintext_email
    outdir
    monochrome_logs
    hook_url
    multiqc_report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: 'nextflow_schema.json')
    def multiqc_reports = multiqc_report.toList()

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error 'Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateInputParameters() {
    genomeExistsError()

    if (params.file == 'fastq-workflow' && !params.input) {
        error("Please provide --input when running the fastq workflow.")
    }
    if (params.file == 'fasta-workflow' && !params.fasta) {
        error("Please provide --fasta when running the fasta workflow.")
    }
    if (params.file == 'primercheck-workflow') {
        if (!params.fasta) {
            error("Please provide --fasta when running the primercheck workflow.")
        }
        if (!params.primer_bed) {
            error("Please provide --primer_bed when running the primercheck workflow.")
        }
        if (!params.primer_fasta) {
            error("Please provide --primer_fasta when running the primercheck workflow.")
        }
    }
}

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}

def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(', ')}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

def toolCitationText() {
    [
        'Tools used in the workflow included:',
        'ARTIC,',
        'FastQC (Andrews 2010),',
        'IRMA (Shepard et al. 2016),',
        'Medaka,',
        'Nextclade (Aksamentov et al. 2021),',
        'and MultiQC (Ewels et al. 2016).',
    ].join(' ').trim()
}

def toolBibliographyText() {
    [
        '<li>Andrews S. FastQC: a quality control tool for high throughput sequence data (2010).</li>',
        '<li>Shepard SS, et al. Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler. BMC Genomics. 2016;17:708.</li>',
        '<li>Aksamentov I, et al. Nextclade: clade assignment, mutation calling and quality control for viral genomes. J Open Source Softw. 2021;6(67):3773.</li>',
        '<li>Ewels P, Magnusson M, Lundin S, Kaller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048.</li>',
    ].join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta['manifest_map'] = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ''
        def manifest_doi = meta.manifest_map.doi.tokenize(',')
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href='https://doi.org/${doi_ref.replace('https://doi.org/', '').replace(' ', '')}'>${doi_ref.replace('https://doi.org/', '').replace(' ', '')}</a>), "
        }
        meta['doi_text'] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else {
        meta['doi_text'] = ''
    }
    meta['nodoi_text'] = meta.manifest_map.doi ? '' : '<li>If available, include the Zenodo DOI for the pipeline release used in the analysis.</li>'
    meta['tool_citations'] = toolCitationText()
    meta['tool_bibliography'] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text
    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

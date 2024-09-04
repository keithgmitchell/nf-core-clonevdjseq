/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_clonevdjseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

help_message = """
    Usage:
        nextflow run nf-core/clonalvdjseq --input <samplesheet.tsv> --rawDataDir <path> --demuxPrimers <demux_primers.fasta> [options]

    Required:
        --input               Samplesheet in TSV format
        --rawDataDir          Directory containing raw sequencing data.
        --demuxPrimers        SMART-Index Barcodes

    Options:
        --outDir              Name of Output directory for results [default: nf-core-clonevdj_results].
        --htstreamOverwrite   Overwrite existing files during HTStream processing [true/false].
        --dada2Overwrite      Overwrite existing files during DADA2 analysis [true/false].
        --help                Print this help message and exit.
    """

// Validate Parameters
if (params.help) {
    println help_message
    exit 0
} else if (!params.input || !params.rawDataDir || !params.demuxPrimers) {
    println help_message
    exit 1, "Please provide all required parameters!"
} 

// BNJ: I hate this solution

if (params.baseDir.endsWith("/")) {
    params.baseDir = params.baseDir[0..-2]
}
if (params.rawDataDir.endsWith("/")) {
    params.rawDataDir = params.rawDataDir[0..-2]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// params for pipeline setup (resources distributed for plates in mount dir)
// get the current directory with pwdnf-core-clonevdj_results using groovy
// params.repobase = System.getProperty("user.dir") + "/"
// params.resourcesDir = params.repobase + "resources/"
// params.samplesheet = params.resourcesDir + "SampleSheet copy.tsv"
// params.metasheet = params.resourcesDir + "SampleSheet copy.tsv"
// params.TSO_demux_primers = params.resourcesDir + "TSO_demux_primers.csv"
// params.baseDir = params.repobase + "exampledata/01-Processing"
// params.hcPrimers = params.resourcesDir + "hc_primers.fasta"
// params.lcPrimers = params.resourcesDir + "lc_primers.fasta"
// params.rawDataDir = params.repobase + "exampledata/00-RawData"

// params for pipeline
params.samplesheet = params.input

// mounting params to use in containers
params.mountDir = params.baseDir + "/" + "nmspipeline"
params.processingDir = params.mountDir + "/" + "01-Processing"
// params.mountedResources = params.mountDir + "resources/"
// params.mtsamplesheet = params.mountedResources + "SampleSheet copy.tsv"
// params.mtmetasheet = params.mountedResources + "alldata_master.tsv"


nextflow.enable.dsl=2


println params.baseDir
println params.samplesheet
println params.baseDir + "/" + params.rawDataDir
println params.baseDir + "/" + params.demuxPrimers


process setupPipeline {
    publishDir "${params.baseDir}/${plate}", mode: 'copy'

    input:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID)

    output:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID)

    script:
    """
    RAW_DATA_DIR=${params.baseDir}/${params.rawDataDir}
    BASE_DIR=${params.baseDir}/${params.outDir}/${plate}

    echo "Setting up pipeline for ${plate}"
    mkdir -p \$BASE_DIR/00-RawData/
    mkdir -p \$BASE_DIR/01-PrimerTrimReport/
    mkdir -p \$BASE_DIR/02-Results/
    mkdir -p \$BASE_DIR/01-PrimerTrim/
    
    r1=\$(find \$RAW_DATA_DIR/ -name '${filePrefix}*_R1_*' | head -n 1)    
    r2=\$(echo \$r1 | sed 's/_R1_/_R2_/')
    echo \$r1
    
    cp \$r1 \$BASE_DIR/00-RawData/
    cp \$r2 \$BASE_DIR/00-RawData/
    #cp ${params.hcPrimers} \$BASE_DIR/hc_primers.fasta
    #cp ${params.lcPrimers} \$BASE_DIR/lc_primers.fasta
    #cp \$NMSEQ_DIR/01-build_hts.py \$BASE_DIR/
    #cp \$NMSEQ_DIR/aberrant_LC.fasta \$BASE_DIR/
    #cp \$NMSEQ_DIR/01-PrimerTrimReport/report.RMD \$BASE_DIR/01-PrimerTrimReport/${plate}_report.RMD
    #cp \$NMSEQ_DIR/SMARTindex_well.tsv \$BASE_DIR/02-Results/
    #cp \$NMSEQ_DIR/02-Results/02-Hybridoma-DADA2-analysis.RMD \$BASE_DIR/02-Results/
    #cp \$NMSEQ_DIR/03-annotate-results.py \$BASE_DIR/
    #cp \$BASE_DIR/03-annotate-results.py .     
    #cp \$NMSEQ_DIR/${Primers} \$BASE_DIR/

    curl ${params.dada2_script} > \$BASE_DIR/02-Results/02-Hybridoma-DADA2-analysis.RMD
    curl ${params.primertrim_script} > \$BASE_DIR/01-PrimerTrimReport/${plate}_report.RMD
    
    #chmod -R 776 ${params.baseDir}
    chmod -R 776 \$BASE_DIR
    """
}

process runHTStream {
    tag "${plate}"
    publishDir "results/${plate}", mode: 'copy'

    input:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID), val(Primer1ID), val(TargetSpecificPrimers), val(TSOBarcode), val(Target_Primer)

    output:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID)

    container 'keithgmitchell/htstream-1.3.3:latest'
    containerOptions '-u $(id -u):$(id -g) --privileged -v ${params.baseDir}:/nmspipeline'

    script:
    """
    BASE_DIR='${params.baseDir}/${params.outDir}/${plate}'
    RAW_DATA_DIR="\${BASE_DIR}/00-RawData"
    LOG_FILE="\${BASE_DIR}/01-PrimerTrim/${TSOBarcode}_${Target_Primer}.log"
    PREFIX="\${BASE_DIR}/01-PrimerTrim/${TSOBarcode}_${Target_Primer}"

    echo "Plate: ${plate}, filePrefix: ${filePrefix}, Primers: ${Primers}, SubmissionID: ${submissionID}, Primer1ID: ${Primer1ID}, TargetSpecificPrimers: ${TargetSpecificPrimers}, TSOBarcode: ${TSOBarcode}, Target_Primer: ${Target_Primer}"
    cd \${BASE_DIR}/00-RawData
    R1_FILES=(\$(ls \${RAW_DATA_DIR}/${filePrefix}*_R1_*))
    R2_FILES=(\$(ls \${RAW_DATA_DIR}/${filePrefix}*_R2_*))
    if [[ \${#R1_FILES[@]} -gt 0 && \${#R2_FILES[@]} -gt 0 ]]; then
        R1=\${R1_FILES[0]}
        R2=\${R2_FILES[0]}
        if [ "${params.htstreamOverwrite}" == "true" ]; then 
            echo "Processing as overwrite is true or no existing files found."
            hts_Primers -d 0 -l 0 -e 0 -r 2 -x -P ${TargetSpecificPrimers} -Q ../${Primer1ID} -1 \${R1} -2 \${R2} -L \${LOG_FILE} |
            hts_NTrimmer -e -A \${LOG_FILE} |
            hts_SeqScreener -C -r -x .01 -k 21 -s ../aberrant_LC.fasta -A \${LOG_FILE} |
            hts_QWindowTrim -l -q 10 -A \${LOG_FILE} |
            hts_Overlapper -A \${LOG_FILE} |
            hts_LengthFilter -m 385 -A \${LOG_FILE} -f \${PREFIX} -F
        else
            echo "HTStream output already exists. Skipping processing."
        fi
    else
        echo "ERROR: Matching R1 or R2 files not found."
    fi
    """
}

process dada2ASVs {
    tag "${plate}"
    publishDir "results/${plate}", mode: 'copy'

    input:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID)
    
    output:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID)

    container 'keithgmitchell/dada2abseq:latest'
    containerOptions '-u $(id -u):$(id -g) --privileged -v ${params.baseDir}:/nmspipeline'

    script:
    """
    BASE_DIR='${params.baseDir}/${params.outDir}/${plate}'
    REPORT_FILE="\${BASE_DIR}/02-Results/02-Hybridoma-DADA2-analysis.html"

    if [ ${params.dada2Overwrite} == "true" ] || [ ! -f "\${REPORT_FILE}" ]; then
        echo "Generating report for ${plate}"
        # read params sheet here in the future?
        Rscript -e "plate='${plate}';mountdir='${params.mountDir}';procdir='${params.processingDir}';submission='${submissionID}';rmarkdown::render('\${BASE_DIR}/02-Results/02-Hybridoma-DADA2-analysis.RMD')"

        cp \$BASE_DIR/02-Results/02-Hybridoma-DADA2-analysis.html ${params.mountDir}/02-Reporting/${plate}_report.html
        # also run the primertrim report
        Rscript -e "plate='${plate}';mountdir='${params.mountDir}';procdir='${params.processingDir}';rmarkdown::render('\${BASE_DIR}/01-PrimerTrimReport/${plate}_report.RMD')"

    else
        echo "DADA2 report already exists. Skipping processing."
    fi
    """
}

process aggregateResults {
    input:
    tuple val(plate), val(filePrefix), val(Primers), val(submissionID)
    
    container 'keithgmitchell/dada2abseq:latest'
    containerOptions '-u $(id -u):$(id -g)'

    script:
    """
    # read params sheet here in the future?
    Rscript -e "metasheet='${params.mtmetasheet}';samplesheet='${params.mtsamplesheet}';resources='${params.mountedResources}';mountdir='${params.mountDir}';procdir='${params.processingDir}';rmarkdown::render('${params.mountedResources}/analyze_plates.rmd')"
    """
}

workflow CLONEVDJSEQ{

    TSOdemux = Channel.fromPath(params.demuxPrimers)
                    .splitCsv(sep: '\t', header: true)
                    .map { row -> tuple(row.TargetSpecificPrimers, row.TSOBarcode, row.Primer1ID, row.Target_Primer) }

    sampleData = Channel.fromPath(params.samplesheet)
                        .splitCsv(sep: '\t', header: true)
                        .map { row -> tuple(row.plate, row.filePrefix, row.Primers, row.submissionID) }

    setupComplete = setupPipeline(sampleData)
    combinedData = setupComplete.combine(TSOdemux)
    htsResults = runHTStream(combinedData)

    groupedByPlate = htsResults.groupTuple(by: [0, 1, 2, 3])
    groupedByPlate.subscribe { groupedData ->
        println "Grouped data: ${groupedData}"
    }

    asvsGenerated = dada2ASVs(groupedByPlate).groupTuple(by: [0, 1, 2, 3]).collect()

    asvsGenerated.subscribe { asvs ->
        println "ASVs generated: ${asvs}"
    }

    aggregateResults(asvsGenerated)

}

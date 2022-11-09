#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/pikavirus
========================================================================================
 nf-core/pikavirus Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/pikavirus
----------------------------------------------------------------------------------------
*/
log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/pikavirus --input samplesheet.csv -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Samplesheet file (-input) not specified!" }

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release']          = workflow.revision
summary['Run Name']                                         = workflow.runName
summary['Input']                                            = params.input
summary['Trimming']                                         = params.trimming
if (params.remove_control) summary['Sequencing control']    = params.control_sequence
summary['Kraken scouting']                                  = params.kraken_scouting
if (params.kraken_scouting) summary['Kraken database']      = params.kraken2_db
summary['Kaiju discovery']                                  = params.kaiju
if (params.translated_analysis) summary ['Kaiju database']  = params.kaiju_db
summary['Virus Search']                                     = params.virus
if (params.virus) summary['Virus Ref']                      = params.vir_ref_dir
if (params.virus) summary['Virus Index File']               = params.vir_dir_repo
summary['Bacteria Search']                                  = params.bacteria
if (params.bacteria) summary['Bacteria Ref']                = params.bact_ref_dir
if (params.bacteria) summary['Bacteria Index File']         = params.bact_dir_repo
summary['Fungi Search']                                     = params.fungi
if (params.fungi) summary['Fungi Ref']                      = params.fungi_ref_dir
if (params.fungi) summary['Fungi Index File']               = params.fungi_dir_repo
summary['Mash winner strategy']                             = params.mash_winner_strategy
summary['Mash identity threshold']                          = params.mash_identity_threshold
summary['Mash min shared hashes']                           = params.mash_min_shared_hashes
summary['Mash pvalue threshold ']                           = params.mash_pvalue_threshold
summary['Max Resources']                                    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container']          = "$workflow.containerEngine - $workflow.container"
summary['Output dir']                                       = params.outdir
summary['Launch dir']                                       = workflow.launchDir
summary['Working dir']                                      = workflow.workDir
summary['Script dir']                                       = workflow.projectDir
summary['User']                                             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']                                   = params.awsregion
    summary['AWS Queue']                                    = params.awsqueue
    summary['AWS CLI']                                      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url

summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']                               = params.email
    summary['E-mail on failure']                            = params.email_on_fail
    summary['MultiQC maxsize']                              = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-pikavirus-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/pikavirus Workflow Summary'
    section_href: 'https://github.com/nf-core/pikavirus'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    errorStrategy 'retry'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }
    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:

    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version &> v_fastqc.txt &
    fastp --version 2> v_fastp.txt
    bowtie2 --version > v_bowtie2.txt
    mash --version &> v_mash.txt &
    spades.py -v &> v_spades.txt &
    quast -v &> v_quast.txt &
    samtools --version | head -n 1 > v_samtools.txt
    bedtools --version > v_bedtools.txt
    kraken2 --version > v_kraken2.txt
    
    kaiju -help &> tmp &
    head -n 1 tmp > v_kaiju.txt

    ivar -v | head -n 1 > v_ivar.txt
    muscle -version > v_muscle.txt

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

process CHECK_SAMPLESHEET {
    tag "$samplesheet"
    publishDir "${params.outdir}/", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.endsWith(".tsv")) "preprocess/sra/$filename"
                      else "pipeline_info/$filename"
    }

    input:
    path(samplesheet) from ch_input

    output:
    path "samplesheet.valid.csv" into ch_samplesheet_reformat
    path "sra_run_info.tsv" optional true

    script:  // These scripts are bundled with the pipeline, in nf-core/viralrecon/bin/
    run_sra = !params.skip_sra && !isOffline()
    """
    awk -F, '{if(\$1 != "" && \$2 != "") {print \$0}}' $samplesheet > nonsra_id.csv
    check_samplesheet.py nonsra_id.csv nonsra.samplesheet.csv
    awk -F, '{if(\$1 != "" && \$2 == "" && \$3 == "") {print \$1}}' $samplesheet > sra_id.list
    if $run_sra && [ -s sra_id.list ]
    then
        fetch_sra_runinfo.py sra_id.list sra_run_info.tsv --platform ILLUMINA --library_layout SINGLE,PAIRED
        sra_runinfo_to_samplesheet.py sra_run_info.tsv sra.samplesheet.csv
    fi
    if [ -f nonsra.samplesheet.csv ]
    then
        head -n 1 nonsra.samplesheet.csv > samplesheet.valid.csv
    else
        head -n 1 sra.samplesheet.csv > samplesheet.valid.csv
    fi
    tail -n +2 -q *sra.samplesheet.csv >> samplesheet.valid.csv
    """
}

// Function to get list of [ sample, single_end?, is_sra?, is_ftp?, [ fastq_1, fastq_2 ], [ md5_1, md5_2] ]
def validate_input(LinkedHashMap sample) {
    def sample_id = sample.sample_id
    def single_end = sample.single_end.toBoolean()
    def is_sra = sample.is_sra.toBoolean()
    def is_ftp = sample.is_ftp.toBoolean()
    def fastq_1 = sample.fastq_1
    def fastq_2 = sample.fastq_2
    def md5_1 = sample.md5_1
    def md5_2 = sample.md5_2

    def array = []
    if (!is_sra) {
        if (single_end) {
            array = [ sample_id, single_end, is_sra, is_ftp, [ file(fastq_1, checkIfExists: true) ] ]
        } else {
            array = [ sample_id, single_end, is_sra, is_ftp, [ file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true) ] ]
        }
    } else {
        array = [ sample_id, single_end, is_sra, is_ftp, [ fastq_1, fastq_2 ], [ md5_1, md5_2 ] ]
    }

    return array
}

/*
 * Create channels for input fastq files
 */
ch_samplesheet_reformat
    .splitCsv(header:true, sep:',')
    .map { validate_input(it) }
    .into { ch_reads_all
            ch_reads_sra }

/*
 * Download and check SRA data
 */
if (!params.skip_sra || !isOffline()) {
    ch_reads_sra
        .filter { it[2] }
        .into { ch_reads_sra_ftp
                ch_reads_sra_dump }

    process SRA_FASTQ_FTP {
        tag "$sample"
        label 'process_medium'
        label 'error_retry'
        publishDir "${params.outdir}/preprocess/sra", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith(".md5")) "md5/$filename"
                          else params.save_sra_fastq ? filename : null
        }

        when:
        is_ftp

        input:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp), val(fastq), val(md5) from ch_reads_sra_ftp

        output:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp), path("*.fastq.gz") into ch_sra_fastq_ftp

        script:
        if (single_end) {
            """
            curl -L ${fastq[0]} -o ${sample}.fastq.gz
            echo "${md5[0]}  ${sample}.fastq.gz" > ${sample}.fastq.gz.md5
            md5sum -c ${sample}.fastq.gz.md5
            """
        } else {
            """
            curl -L ${fastq[0]} -o ${sample}_1.fastq.gz
            echo "${md5[0]}  ${sample}_1.fastq.gz" > ${sample}_1.fastq.gz.md5
            md5sum -c ${sample}_1.fastq.gz.md5
            curl -L ${fastq[1]} -o ${sample}_2.fastq.gz
            echo "${md5[1]}  ${sample}_2.fastq.gz" > ${sample}_2.fastq.gz.md5
            md5sum -c ${sample}_2.fastq.gz.md5
            """
        }
    }

    process SRA_FASTQ_DUMP {
        tag "$sample"
        label 'process_medium'
        label 'error_retry'
        publishDir "${params.outdir}/preprocess/sra", mode: params.publish_dir_mode,
            saveAs: { filename ->
                          if (filename.endsWith(".log")) "log/$filename"
                          else params.save_sra_fastq ? filename : null
        }

        when:
        !is_ftp

        input:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp) from ch_reads_sra_dump.map { it[0..3] }

        output:
        tuple val(sample), val(single_end), val(is_sra), val(is_ftp), path("*.fastq.gz") into ch_sra_fastq_dump
        path "*.log"

        script:
        prefix = "${sample.split('_')[0..-2].join('_')}"
        pe = single_end ? "" : "--readids --split-e"
        rm_orphan = single_end ? "" : "[ -f  ${prefix}.fastq.gz ] && rm ${prefix}.fastq.gz"
        """
        parallel-fastq-dump \\
            --sra-id $prefix \\
            --threads $task.cpus \\
            --outdir ./ \\
            --tmpdir ./ \\
            --gzip \\
            $pe \\
            > ${prefix}.fastq_dump.log
        $rm_orphan
        """
    }

    ch_reads_all
        .filter { !it[2] }
        .concat(ch_sra_fastq_ftp, ch_sra_fastq_dump)
        .set { ch_reads_all }
}

ch_reads_all
    .map { [ it[0].split('_')[0..-2].join('_'), it[1], it[4] ] }
    .groupTuple(by: [0, 1])
    .map { [ it[0], it[1], it[2].flatten() ] }
    .set { ch_reads_all }


/*
 * Merge FastQ files with the same sample identifier (resequenced samples)
 */

process CAT_FASTQ {
    tag "$sample"

    input:
    tuple val(sample), val(single_end), path(reads) from ch_reads_all

    output:
    tuple val(sample), val(single_end), path("*.merged.fastq.gz") into ch_cat_fastqc,
                                                                       ch_cat_fortrim
    val(sample) into samplechannel_translated,
                         samplechannel_trim_fastqc,
                         samplechannel_trimmed_fastqc,
                         samplechannel_krona,
                         control_results_template,
                         virus_results_template,
                         bacteria_results_template,
                         fungi_results_template

    tuple val(sample), val(single_end) into sample_pe_template

    script:
    readList = reads.collect{it.toString()}
    if (!single_end) {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.sort().join(' ')} > ${sample}_1.merged.fastq.gz
            cat ${read2.sort().join(' ')} > ${sample}_2.merged.fastq.gz
            """
        } else {
            """
            ln -s ${reads[0]} ${sample}_1.merged.fastq.gz
            ln -s ${reads[1]} ${sample}_2.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 1) {
            """​​​​​​​
            cat ${readList.sort().join(' ')} > ${sample}.merged.fastq.gz
            """
        } else {
            """
            ln -s $reads ${sample}.merged.fastq.gz
            """
        }
    }
}


/*
 * PREPROCESSING: KAIJU DATABASE
 */
if (params.kaiju && params.translated_analysis) {

    if (params.kaiju_db.endsWith('.gz') || params.kaiju_db.endsWith('.tar') || params.kaiju_db.endsWith('.tgz')){

        process UNCOMPRESS_KAIJUDB {
            label 'error_retry'

            input:
            path(database) from params.kaiju_db

            output:
            path("kaijudb") into kaiju_db

            script:
            """
            mkdir "kaijudb"
            tar -zxf $database -C "kaijudb"
            """
        }
    

    } else {
        kaiju_db = Channel.fromPath(params.kaiju_db)
    }
}

/*
 * STEP 1.1 - FastQC
 */
process RAW_SAMPLES_FASTQC {
    tag "$samplename"
    label "process_medium"
    publishDir "${params.outdir}/${samplename}/raw_fastqc", mode: params.publish_dir_mode,
    saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
    }

    input:
    set val(samplename), val(single_end), path(reads) from ch_cat_fastqc

    output:
    tuple val(samplename), path("*_fastqc.{zip,html}")

    tuple val(samplename), path("*_fastqc.zip") into raw_fastqc_multiqc
    path("*_fastqc.zip") into raw_fastqc_multiqc_global

    script:

    """
    fastqc --quiet --threads $task.cpus $reads
    """
}

/*
 * STEP 1.2 - TRIMMING
*/
if (params.trimming) {
    process FASTP_TRIM {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}", mode: params.publish_dir_mode,
            saveAs: { filename ->
                        if (filename.endsWith(".fastq") && params.rescue_trimmed) "trimmed_sequences/$filename"
                        else if (filename.endsWith(".html")) filename
                    }

        input:
        tuple val(samplename), val(single_end), path(reads) from ch_cat_fortrim

        output:
        tuple val(samplename), val(single_end), path("*fail.fastq.gz")

        tuple val(samplename), val(single_end), path("*trim.fastq.gz") into trimmed_paired_fastqc, trimmed_remove_control
        tuple val(samplename), path("*.json") into fastp_multiqc
        tuple val(samplename), path("*.html") into fastp_report
        path("*.json") into fastp_multiqc_global

        script:
        detect_adapter =  single_end ? "" : "--detect_adapter_for_pe"
        reads1 = single_end ? "--in1 ${reads} --out1 ${samplename}_trim.fastq.gz --failed_out ${samplename}_fail.fastq.gz" : "--in1 ${reads[0]} --out1 ${samplename}_1_trim.fastq.gz --unpaired1 ${samplename}_1_fail.fastq.gz"
        reads2 = single_end ? "" : "--in2 ${reads[1]} --out2 ${samplename}_2_trim.fastq.gz --unpaired2 ${samplename}_2_fail.fastq.gz"

        """
        fastp \\
        $detect_adapter \\
        --cut_front \\
        --cut_tail \\
        --length_required 35 \\
        --thread $task.cpus \\
        --json ${samplename}_trim_fastp.json \\
        $reads1 \\
        $reads2
        """
    }

    /*
    * STEP 1.3 - FastQC on trimmed reads
    */
    process TRIMMED_SAMPLES_FASTQC {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/trimmed_fastqc", mode: params.publish_dir_mode

        input:
        tuple val(samplename), val(single_end), path(reads) from trimmed_paired_fastqc

        output:
        tuple val(samplename), path("*_fastqc.{zip,html}")
        tuple val(samplename), path("*_fastqc.zip") into trimmed_fastqc_multiqc
        path("*_fastqc.zip") into trimmed_fastqc_multiqc_global

        script:

        """
        fastqc --quiet --threads $task.cpus $reads
        """
    }

    if (params.remove_control) {

        Channel.fromPath(params.control_sequence).set { control_genome }

        process BOWTIE2_REMOVE_SEQUENCING_CONTROL {
            tag "$samplename"
            label "process_high"

            input:
            tuple val(samplename), val(single_end), path(reads), path(control_sequence) from trimmed_remove_control.combine(control_genome)
            
            output:
            tuple val(samplename), val(single_end), path("*_mapped_sorted.bam") into control_alignment
            tuple val(samplename), val(single_end), path("*.fastq.gz") into trimmed_virus, trimmed_bact, trimmed_fungi,
                                                                            reads_for_assembly, trimmed_skip_hostremoval,
                                                                            trimmed_map_virus, trimmed_map_bact, trimmed_map_fungi,
                                                                            trimmed_kraken2
            script:
            samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
            unmapped = single_end ? "--un-gz ${samplename}_unmapped.fastq.gz" : "--un-conc-gz ${samplename}_unmapped_R%.fastq.gz"

            """
            bowtie2-build \\
            --seed 1 \\
            --threads $task.cpus \\
            $control_sequence \\
            "control_sequence"

            bowtie2 \\
            --threads $task.cpus \\
            -x "control_sequence" \\
            $unmapped \\
            $samplereads | samtools view \\
            -@ $task.cpus \\
            -b \\
            -h \\
            -O BAM \\
            -o "${control_sequence}_mapped.bam"

            samtools sort \\
            -@ $task.cpus \\
            -o "${control_sequence}_mapped_sorted.bam" \\
            "${control_sequence}_mapped.bam"

            """
        }

        process SAMTOOLS_CONTROL {
            tag "$samplename"
            label "process_medium"

            input:
            tuple val(samplename), val(single_end), path(sortedbam) from control_alignment

            output:
            tuple val(samplename), val(single_end), path(sortedbam), path("*idxstats"), path("*flagstat") into control_alignment_bams

            script:

            prefix = sortedbam.join().minus("_mapped_sorted.bam")

            """
            samtools index $sortedbam
            samtools idxstats $sortedbam > "${prefix}.sorted.bam.idxstats"
            samtools flagstat -O tsv $sortedbam > "${prefix}.sorted.bam.flagstat"

            """
        }

        process BEDTOOLS_SEQUENCING_CONTROL {
            tag "$samplename"
            label "process_medium"


            input:
            tuple val(samplename), val(single_end), path(mapped), path(idxstat), path(flagstat) from control_alignment_bams

            output:
            tuple val(samplename), path("*_coverage.txt"), path(idxstat), path(flagstat) into coverage_stats_control

            script:
            prefix = mapped.join().minus("_mapped_sorted.bam")

            """
            bedtools genomecov -ibam $mapped  > "${prefix}_coverage.txt"
            """

        }

        process COVERAGE_STATS_SEQUENCING_CONTROL {
            tag "$samplename"
            label "process_medium"

            input:
            tuple val(samplename), path(coveragefile), path(idxstat), path(flagstat) from coverage_stats_control

            output:
            tuple val(samplename), path("*.tsv") into control_coverage_results

            script:
            """
            coverage_analysis_control.py $samplename $coveragefile $idxstat $flagstat
            """
        }

    } else {
        trimmed_remove_control.into { trimmed_virus
                                      trimmed_bact
                                      trimmed_fungi
                                      reads_for_assembly
                                      trimmed_skip_hostremoval
                                      trimmed_map_virus
                                      trimmed_map_bact
                                      trimmed_map_fungi
                                      trimmed_kraken2 }

        nofile_path_control_coverage = Channel.fromPath("NONE_control_coverage")
        control_results_template.combine(nofile_path_control_coverage).set { control_coverage_results }

    }

} else {

    ch_cat_fortrim.into { trimmed_virus
                          trimmed_bact
                          trimmed_fungi
                          reads_for_assembly
                          trimmed_skip_hostremoval
                          trimmed_map_virus
                          trimmed_map_bact
                          trimmed_map_fungi }

    trimmed_fastqc_multiqc_global = Channel.fromPath("NONE_trimmedfastqc_global")
    fastp_multiqc_global = Channel.fromPath("NONE_fastp_global")

    nofile_path_trimmed_fastqc = Channel.fromPath("NONE_trimmedfastqc")
    nofile_path_fastp = Channel.fromPath("NONE_fastp")

    nofile_path_control_coverage = Channel.fromPath("NONE_control_coverage")
    control_results_template.combine(nofile_path_control_coverage).set { control_coverage_results }

    samplechannel_trim_fastqc.combine(nofile_path_fastp).set { fastp_multiqc }
    samplechannel_trimmed_fastqc.combine(nofile_path_trimmed_fastqc).set { trimmed_fastqc_multiqc }

}

if (params.kraken_scouting || params.translated_analysis) {

    if (params.kraken2_db.endsWith('.gz') || params.kraken2_db.endsWith('.tar') || params.kraken2_db.endsWith('.tgz')) {

        Channel.fromPath(params.kraken2_db).set { kraken2_compressed }

        process UNCOMPRESS_KRAKEN2DB {
            label 'error_retry'

            input:
            path(database) from kraken2_compressed

            output:
            path("kraken2db") into kraken2_db_files

            script:
            """
            tar -zxf $database

            if [ -f "hash.k2d" ];
            then
                mkdir kraken2db

                mv *.k2d kraken2db
                mv *.kmer_distrib kraken2db

                mv seqid2taxid.map kraken2db
                mv inspect.txt kraken2db
                mv README_assembly_summary.txt kraken2db

            else
                mv */ kraken2db
            fi
            
            """
        }
    } else {
        Channel.fromPath(params.kraken2_db).set { kraken2_db_files }
    }

    process SCOUT_KRAKEN2 {
        tag "$samplename"
        label "process_high"

        input:
        tuple val(samplename), val(single_end), path(reads), path(kraken2db) from trimmed_kraken2.combine(kraken2_db_files)

        output:
        tuple val(samplename), path("*.report") 
        tuple val(samplename), path("*.krona") into kraken2_krona
        tuple val(samplename), path("*.report"), path("*.kraken") into kraken2_host
        tuple val(samplename), val(single_end), file("*_unclassified.fastq") into unclassified_reads

        script:
        paired_end = single_end ? "" : "--paired"
        unclass_name = single_end ? "${samplename}_unclassified.fastq" : "${samplename}_#_unclassified.fastq"
        """
        kraken2 --db $kraken2db \\
        ${paired_end} \\
        --threads $task.cpus \\
        --report ${samplename}.report \\
        --output ${samplename}.kraken \\
        --unclassified-out ${unclass_name} \\
        ${reads}
        cat ${samplename}.kraken | cut -f 2,3 > results.krona
        """
    }

    if (params.kraken_scouting) {

        
        process KRONA_DB {

            output:
            path("taxonomy/") into krona_taxonomy_db_kraken

            script:
            """
            ktUpdateTaxonomy.sh taxonomy
            """
        }

        process KRONA_KRAKEN_RESULTS {
            tag "$samplename"
            label "process_medium"
            publishDir "${params.outdir}/${samplename}/kraken2_krona_results", mode: params.publish_dir_mode

            input:
            tuple val(samplename), path(kronafile), path(taxonomy) from kraken2_krona.combine(krona_taxonomy_db_kraken)

            output:
            tuple val(samplename), path("*.krona.html") into krona_scouting_results

            script:
            outfile = "${samplename}_kraken.krona.html"
            
            """
            ktImportTaxonomy $kronafile -tax $taxonomy -o $outfile
            """
        }

    }

    if (params.host_removal) {

        process REMOVE_HOST_KRAKEN2 {
            tag "$samplename"
            label "process_medium"

            input:
            tuple val(samplename), val(single_end), path(reads), path(report), path(output) from trimmed_paired_extract_host.join(kraken2_host_extraction)

            output:
            tuple val(samplename), val(single_end), path("*_host_extracted.fastq") into reads_for_assembly

            script:
            read = single_end ? "-s ${reads}" : "-s1 ${reads[0]} -s2 ${reads[1]}"
            outputfile = single_end ? "--output $mergedfile" : "-o ${samplename}_1_host_extracted.fastq -o2 ${samplename}_2_host_extracted.fastq"
            mergedfile = single_end ? "${samplename}_host_extracted.fastq": "${samplename}_merged.fastq"
            merge_outputfile = single_end ? "" : "cat *_host_extracted.fastq > $mergedfile"
            """
            extract_kraken_reads.py \\
            -k $output \\
            -r $report \\
            --exclude \\
            --taxid ${params.host_taxid} \\
            --fastq-output \\
            $read \\
            $outputfile
            $merge_outputfile
            """
        }

    } else {
        trimmed_skip_hostremoval.set { reads_for_assembly }
    }

} else {

    Channel.fromPath("NONE_krona").set { nofile_path_krona }
    samplechannel_krona.combine(nofile_path_krona).set { krona_scouting_results }

}

if (params.virus) {

    Channel.fromPath(params.vir_dir_repo).into { virus_datasheet_coverage
                                                 virus_datasheet_len
                                                 virus_datasheet_selection
                                                 virus_datasheet_group_by_species }

    if (params.vir_ref_dir.endsWith('.gz') || params.vir_ref_dir.endsWith('.tar') || params.vir_ref_dir.endsWith('.tgz')) {

        process UNCOMPRESS_VIRUS_REF {
            label 'error_retry'

            input:
            path(ref_vir) from params.vir_ref_dir

            output:
            path("viralrefs") into virus_ref_directory, virus_references

            script:
            """
            mkdir "viralrefs"
            tar -xvf $ref_vir --strip-components=1 -C "viralrefs"
            """
        }

    } else {
        Channel.fromPath(params.vir_ref_dir).into { virus_ref_directory
                                                    virus_references }
    }

    process MASH_GENERATE_REFERENCE_SKETCH_VIRUS {
        label "process_high"

        input:
        path(ref) from virus_ref_directory

        output:
        path("*.msh") into reference_sketch_virus

        script:
        """
        find ${ref}/ -name "*" -type f > reference_list.txt
        mash sketch -k 32 -s 5000 -o reference -l reference_list.txt
        """

    }

    process MASH_DETECT_VIRUS_REFERENCES {
        tag "$samplename"
        label "process_high"
        publishDir "${params.outdir}/mash_results", mode: params.publish_dir_mode

        input:
        tuple val(samplename), val(single_end), path(reads), path(refsketch) from trimmed_virus.combine(reference_sketch_virus)

        output:
        tuple val(samplename), path(mashout) into mash_result_virus_references

        script:
        mashout = "mash_screen_results_virus_${samplename}.txt"
        winstrat = params.mash_winner_strategy ? "-w" : ""
        
        """
        echo -e "#Identity\tShared_hashes\tMedian_multiplicity\tP-value\tQuery_id\tQuery_comment" > $mashout
        mash screen $winstrat $refsketch $reads >> $mashout
        """
    }

    process SELECT_FINAL_VIRUS_REFERENCES {
        tag "$samplename"
        label "process_low"
        publishDir "${params.outdir}/${samplename}/virus_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.endsWith(".tsv")) filename
        }

        input:
        tuple val(samplename), path(mashresult), path(refdir), path(datasheet_virus) from mash_result_virus_references.combine(virus_references).combine(virus_datasheet_selection)

        output:
        tuple val(samplename), path("Final_fnas/*") optional true into bowtie_virus_references
        tuple val(samplename), path("not_found.tsv") optional true into failed_virus_samples

        script:
        """
        extract_significative_references.py \\
        --mash-result $mashresult \\
        --refdir $refdir \\
        --ref-sheet $datasheet_virus \\
        --identity-threshold $params.mash_identity_threshold \\
        --shared-hashes-threshold $params.mash_min_shared_hashes \\
        --p-value-threshold $params.mash_pvalue_threshold
        """
    }

    
    trimmed_map_virus.join(bowtie_virus_references).set { bowtie_virus_channel }

    // Channel is: [ val(samplesheet), path(reads), path(a lot of files) ]
    // the following code turns it into [ val(samplesheet), path (reads), path(only one single file)]

    def rawlist_virus = bowtie_virus_channel.toList().get()
    def bowtielist_virus = []

    for (line in rawlist_virus) {
        if (line[3] instanceof java.util.ArrayList){
            last_list = line[3]
            }
            else {
                last_list = [line[3]]
            }

            for (reference in last_list) {
                def ref_slice = [line[0],line[1],line[2],reference]
                bowtielist_virus.add(ref_slice)
        }
    }

    def reads_virus_mapping = Channel.fromList(bowtielist_virus)

    process BOWTIE2_MAPPING_VIRUS {
        tag "${samplename} : ${reference_sequence}"
        label "process_high"

        input:
        tuple val(samplename), val(single_end), path(reads), path(reference_sequence) from reads_virus_mapping

        output:
        tuple val(samplename), val(single_end), path("*sorted.bam") into bowtie_alingment_bam_virus, bowtie_alingment_bam_virus_ivar
        tuple val(samplename), val(single_end), path("*sorted.bam"), path(reference_sequence) into ivar_virus
        tuple val(samplename), path("*.sam") into bowtie_alignment_sam_virus

        script:
        samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
        prefix = "${reference_sequence}_vs_${samplename}_virus"
        
        """
        bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $reference_sequence \\
        "index"

        bowtie2 \\
        --threads $task.cpus \\
        -x "index" \\
        $samplereads > ${prefix}.sam

        samtools view \\
        -@ $task.cpus \\
        -b \\
        -h \\
        -O BAM \\
        -o "${prefix}.bam" \\
        "${prefix}.sam"

        samtools sort \\
        -@ $task.cpus \\
        -o "${prefix}.sorted.bam" \\
        "${prefix}.bam"

        """
    }

    process SAMTOOLS_VIRUS {
        tag "${samplename} : ${reference_sequence}"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(sortedbam) from bowtie_alingment_bam_virus_ivar

        output:
        tuple val(samplename), val(single_end), path("*.flagstat"), path("*.idxstats"), path("*.stats")
        
        script:
        prefix = sortedbam.join()

        """
        samtools index $sortedbam

        samtools flagstat -O tsv $sortedbam > "${prefix}.flagstat"
        samtools idxstats $sortedbam > "${prefix}.idxstats"
        samtools stats $sortedbam > "${prefix}.stats"
        """
    }

    process IVAR_CONSENSUS_VIRUS {
        tag "${samplename} : ${prefix}"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/virus_consensus"
        
        input:
        tuple val(samplename), val(single_end), path(sortedbam), path(reference_sequence) from ivar_virus

        output:
        tuple val(samplename), path("*.fa") into ch_ivar_consensus

        script:
        prefix = sortedbam.join().minus("_virus.sorted.bam")

        """

        if [[ $reference_sequence == **.gz ]]
        then
            gunzip -c $reference_sequence > fastaref

        else

            mv $reference_sequence fastaref
        fi

        samtools mpileup \\
            --count-orphans \\
            --no-BAQ \\
            --fasta-ref fastaref \\
            --min-BQ 20 \\
            --output ${prefix}.mpileup \\
            $sortedbam

        rm -rf fastaref

        cat ${prefix}.mpileup | ivar consensus -t 0.51 -n N -p ${prefix}_consensus

        """
    }

    process GROUP_BY_ORGANISM_VIRUS {
        tag "${samplename}"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/consensus_sequences_virus", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.contains("_consensus_sequence")) filename
                    }

        input:
        tuple val(samplename), path(consensus_files), path(datasheet_virus) from ch_ivar_consensus.groupTuple().combine(virus_datasheet_group_by_species)

        output:
        tuple val(samplename), path("*_consensus_sequence*") optional true
        tuple val(samplename), path("*_directory") optional true into virus_consensus_by_species_raw

        script:
        """
        organism_attribution.py $samplename $datasheet_virus $consensus_files
        """

    }

    def virus_species_consensus_list = virus_consensus_by_species_raw.toList().get()
    if (virus_species_consensus_list.size() > 0) {

        def consensus_list = []

        for (sample in virus_species_consensus_list) {

            def samplename = sample[0]

            if (sample[1] instanceof java.util.ArrayList) {

                for (consensus_dir in sample[1]) {
                    def consensus_slice = [samplename, consensus_dir]
                    consensus_list.add(consensus_slice)
                }

            } else {

                def consensus_slice = [samplename, sample[1]]
                consensus_list.add(consensus_slice)

            }
        }

        Channel.fromList(consensus_list).set { virus_consensus_by_species }

    } else {

        Channel.empty().set { virus_consensus_by_species }
    }



    process MUSCLE_ALIGN_CONSENSUS_VIRUS {
        tag "${samplename}: ${prefix}"
        label "process_high"
        publishDir "${params.outdir}/${samplename}/consensus/", mode: params.publish_dir_mode


        input:
        tuple val(samplename), path(consensus_dir) from virus_consensus_by_species

        output:
        tuple val(samplename), path("*_msa.fasta")

        script:
        prefix = consensus_dir.join().minus("_consensus_directory").minus("/")
        """

        cat ${consensus_dir}/* > multifasta

        muscle -in multifasta \\
               -maxiters 2 \\ 
               -out ${prefix}_msa.fasta

        """
    }

    process BEDTOOLS_COVERAGE_VIRUS {
        tag "$samplename"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(bamfiles) from bowtie_alingment_bam_virus

        output:
        tuple val(samplename), path("*_bedgraph_virus.txt") into bedgraph_virus
        tuple val(samplename), path("*_coverage_virus.txt") into coverage_files_virus_merge

        script:
        prefix = bamfiles.join().minus("sorted.bam")
        """
        bedtools genomecov -ibam $bamfiles  > "${prefix}_coverage_virus.txt"
        bedtools genomecov -ibam $bamfiles -bga >"${prefix}_bedgraph_virus.txt"
        """
    }

    process COVERAGE_STATS_VIRUS {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.endsWith(".html")) "virus_coverage/plots/$filename"
                      else if (filename.endsWith(".tsv")) filename
                      else "virus_coverage/$filename"
                    }


        input:
        tuple val(samplename), path(coveragefiles), path(datasheet_virus) from coverage_files_virus_merge.groupTuple().combine(virus_datasheet_coverage)

        output:
        path("*.html")
        tuple val(samplename), path("*.tsv") into coverage_stats_virus
        path("*.tsv") into coverage_stats_tomerge_virus
        path("*_valid_coverage_files_virus") into valid_coverage_files_virus

        script:
        """
        graphs_coverage.py $samplename virus $datasheet_virus $coveragefiles
        """
    }


    // Channel is: [ val(samplename), path("a single sam file") ]
    // Turn it into: [ val(samplename), path("all sam files in for the sample")]

    def tmp_list  = bowtie_alignment_sam_virus.toList().get()

    def sam_virus = [:]
    
    for (line in tmp_list) {
        if (sam_virus.containsKey(line[0])) {
            sam_virus[line[0]].add(line[1])
        } else {
            sam_virus[line[0]] = [line[1]]
        }
    }
    
    def organized_sam_list_virus = []
    for (entry in sam_virus) {
        slice = [entry.key]
        slice.add(entry.value)
        organized_sam_list_virus.add(slice)
    }

    def ch_sam_virus = Channel.fromList(organized_sam_list_virus)

    process FIND_UNIQUE_READS_VIRUS {
        tag "$samplename"
        label "process_low"

        input:
        tuple val(samplename), path(samfiles) from ch_sam_virus

        output:
        tuple val(samplename), path(samfiles), path("*_mapping_balance.tsv")
        tuple val(samplename), path(samfiles), path("*_mapped_reads.txt") into mapped_reads
        tuple val(samplename), path(samfiles), path("*_unmapped_reads.txt") into unmapped_reads
        tuple val(samplename), path(samfiles), path("*_unique_reads.txt") into unique_reads
        
        script:
        """
        find_unique_reads.py ${samplename}
        """
    }


    if (params.keep_mapped_reads_bam == true) {
    // Channel is  : [ val(samplename), path("a lot of sam files"), path("a lot of txt files")]
    // Turn it into: [ val(samplename), path("a single sam file" ), path("a single txt file" )] 
        
        def mapped_list       = mapped_reads.toList().get()
        def mapped_reads_list = []

        for ( item in mapped_list ) {
            def mapped_map = [:]
            
            if ( item[1] instanceof java.util.ArrayList ) {
                sam_list = item[1]
            } else {
                sam_list = [item[1]]
            }

            for ( sam_file in sam_list ) {
                assembly_name = sam_file.last().toString() - item[0] - "_vs_" - ".sam" -".fna.gz" - "_virus" - "_bacteria" - "_fungi"
                if (! mapped_map.containsKey(assembly_name)) {
                    mapped_map[assembly_name] = [sam_file]
                }
            }

            if ( item[2] instanceof java.util.ArrayList ) {
                report_list = item[2]
            } else {
                report_list = [item[2]]
            }

            for ( mapped_report in report_list ) {
                for ( key in mapped_map.keySet()) {
                    if ( mapped_report.toString().contains(key) ) {
                        mapped_map[key].add(mapped_report)
                    }
                }
            }

            for (entry in mapped_map) {
                slice = [item[0], entry.key, entry.value[0], entry.value[1]]
                mapped_reads_list.add(slice)
            }
        }
        
        def mapped_reads = Channel.fromList(mapped_reads_list)

        process EXTRACT_MAPPED_READS {
            tag "$samplename"
            label "process_medium"
            publishDir "${params.outdir}/${samplename}/", mode: params.publish_dir_mode

            input:
            tuple val(samplename), val(assembly_name), path(samfile), path(report_list) from mapped_reads

            script:
            filtered_sam = "mapped_reads_${assembly_name}.sam"
            """
            grep -e "^@" ${samfile} > ${filtered_sam}
            cat ${report_list} | xargs -I @@ grep @@ ${samfile} >> ${filtered_sam}
            """
    
        }



    }

    // if (params.keep_unique_reads_bam == true) {}

    // if (params.keep_unmapped_reads_bam == true) {}

    process MERGE_COVERAGE_TABLES_VIRUS {
        label "process_low"
        publishDir "${params.outdir}", mode: params.publish_dir_mode

        input:
        path(coverage_tsvs) from coverage_stats_tomerge_virus.collect()

        output:
        path("all_samples_virus_table.tsv")

        script:
        """
        echo -e "samplename\tgnm\tspecies\tsubspecies\tcovMean\tcovSD\tcovMin\tcovMax\tcovMedian\t>=x1\t>=x10\t>=x25\t>=x50\t>=x75\t>=x100\tassembly" > all_samples_virus_table

        for item in ./*.tsv;
        do
            samplename=\$(echo \$item | sed s/"_virus_table.tsv"// | sed s@"./"@@)
            awk -F "\t" -v name="\$samplename" 'BEGIN { OFS = FS } { if (NR>1) { \$1=name; print } }' \$item >> all_samples_virus_table
        done

        mv all_samples_virus_table all_samples_virus_table.tsv

        """
    }

    process COVERAGE_LEN_VIRUS {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/virus_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.endsWith(".html")) "plots/$filename"
                      else filename
        }

        input:
        tuple val(samplename), path(bedgraph), path(datasheet_virus) from bedgraph_virus.groupTuple().combine(virus_datasheet_len)

        output:
        path("*.html")
        path("*_valid_bedgraph_files_virus") into valid_bedgraph_files_virus

        script:
        """
        generate_len_coverage_graph.py $samplename virus $datasheet_virus $bedgraph
        """
    }

    coverage_stats_virus.mix(failed_virus_samples).set { virus_coverage_results }

} else {

    nofile_path_virus_coverage = Channel.fromPath("NONE_virus")
    virus_results_template.combine(nofile_path_virus_coverage).set { virus_coverage_resul }

}

if (params.bacteria) {

    Channel.fromPath(params.bact_dir_repo).into { bact_table
                                                  bact_table_len
                                                  bact_sheet }

    if (params.bact_ref_dir.endsWith('.gz') || params.bact_ref_dir.endsWith('.tar') || params.bact_ref_dir.endsWith('.tgz')) {

        process UNCOMPRESS_BACTERIA_REF {
            label 'error_retry'

            input:
            path(ref_bact) from params.bact_ref_dir

            output:
            path("bactrefs") into bact_ref_directory, bact_references

            script:
            """
            mkdir "bactrefs"
            tar -xvf $ref_bact --strip-components=1 -C "bactrefs"
            """
        }
    } else {
        Channel.fromPath(params.bact_ref_dir).into { bact_ref_directory
                                                     bact_references }
    }

    process MASH_DETECT_BACTERIA_REFERENCES {
        tag "$samplename"
        label "process_high"

        input:
        tuple val(samplename), val(single_end), path(reads), path(ref) from trimmed_bact.combine(bact_ref_directory)

        output:
        tuple val(samplename), path(mashout) into mash_result_bact_references

        script:
        mashout = "mash_screen_results_bact_${samplename}.txt"

        """
        find ${ref}/ -name "*" -type f > reference_list.txt
        mash sketch -k 32 -s 5000 -o reference -l reference_list.txt
        echo -e "#Identity\tShared_hashes\tMedian_multiplicity\tP-value\tQuery_id\tQuery_comment" > $mashout
        mash screen reference.msh $reads >> $mashout
        """
    }

    process SELECT_FINAL_BACTERIA_REFERENCES {
        tag "$samplename"
        label "process_low"
        publishDir "${params.outdir}/${samplename}/bacteria_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.endsWith(".tsv")) filename
        }

        input:
        tuple val(samplename), path(mashresult), path(refdir), path(datasheet) from mash_result_bact_references.combine(bact_references).combine(bact_sheet)

        output:
        tuple val(samplename), path("Final_fnas/*") optional true into bowtie_bact_references
        tuple val(samplename), path("not_found.tsv") optional true into failed_bact_samples

        script:
        """
        extract_significative_references.py $mashresult $refdir $datasheet

        """
    }

    trimmed_map_bact.join(bowtie_bact_references).set { bowtie_bact_channel }

    def rawlist_bact = bowtie_bact_channel.toList().get()
    def bowtielist_bact = []

    for (line in rawlist_bact) {
        if (line[3] instanceof java.util.ArrayList){
            last_list = line[3]
            }
            else {
                last_list = [line[3]]
            }

            for (reference in last_list) {
                def ref_slice = [line[0],line[1],line[2],reference]
                bowtielist_bact.add(ref_slice)
        }
    }

    def reads_bact_mapping = Channel.fromList(bowtielist_bact)

    process BOWTIE2_MAPPING_BACTERIA {
        tag "${samplename} : ${reference}"
        label "process_high"

        input:
        tuple val(samplename), val(single_end), path(reads), path(reference) from reads_bact_mapping

        output:
        tuple val(samplename), val(single_end), path("*_bact.sam"), path(reference) into bowtie_alingment_sam_bact

        script:
        samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
        """
        bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $reference \\
        "index"

        bowtie2 \\
        -x "index" \\
        ${samplereads} \\
        -S "${reference}_vs_${samplename}_bact.sam" \\
        --threads $task.cpus
        """
    }

   process SAMTOOLS_BACTERIA {
        tag "$samplename"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(samfiles), path(reference) from bowtie_alingment_sam_bact

        output:
        tuple val(samplename), val(single_end), path("*.sorted.bam") into bowtie_alingment_bam_bact
        tuple val(samplename), val(single_end), path("*.sorted.bam"), path(reference) into ordered_bam_mpileup_bact
        tuple val(samplename), val(single_end), path("*.sorted.bam.flagstat"), path("*.sorted.bam.idxstats"), path("*.sorted.bam.stats") into bam_stats_bact

        script:

        """
        samtools view \\
        -@ $task.cpus \\
        -b \\
        -h \\
        -O BAM \\
        -o "\$(basename $samfiles .sam).bam" \\
        $samfiles

        samtools sort \\
        -@ $task.cpus \\
        -o "\$(basename $samfiles .sam).sorted.bam" \\
        "\$(basename $samfiles .sam).bam"
        samtools index "\$(basename $samfiles .sam).sorted.bam"
        samtools flagstat "\$(basename $samfiles .sam).sorted.bam" > "\$(basename $samfiles .sam).sorted.bam.flagstat"
        samtools idxstats "\$(basename $samfiles .sam).sorted.bam" > "\$(basename $samfiles .sam).sorted.bam.idxstats"
        samtools stats "\$(basename $samfiles .sam).sorted.bam" > "\$(basename $samfiles .sam).sorted.bam.stats"
        """
    }

    process BEDTOOLS_COVERAGE_BACTERIA {
        tag "$samplename"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(bamfiles) from bowtie_alingment_bam_bact

        output:
        tuple val(samplename), path("*_bedgraph.txt") into bedgraph_bact
        tuple val(samplename), path("*_coverage.txt") into coverage_files_bact_merge


        script:

        """
        bedtools genomecov -ibam $bamfiles  > "\$(basename -- $bamfiles .sorted.bam)_coverage.txt"
        bedtools genomecov -ibam $bamfiles  -bga >"\$(basename -- $bamfiles .sorted.bam)_bedgraph.txt"
        """
    }

    process COVERAGE_STATS_BACTERIA {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/bacteria_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                     if (filename.endsWith(".html")) "plots/$filename"
                     else "$filename"

        }

        input:
        tuple val(samplename), path(coveragefiles), path(reference_bacteria) from coverage_files_bact_merge.groupTuple().combine(bact_table)

        output:
        tuple val(samplename), path("*.tsv") into coverage_stats_bacteria
        path("*.html") into coverage_graphs_bacteria
        path("*_valid_coverage_files_virus") into valid_coverage_files_bacteria

        script:

        """
        graphs_coverage.py $samplename bacteria $reference_bacteria $coveragefiles
        """
    }

    process COVERAGE_LEN_BACTERIA {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/bacteria_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.endsWith(".html")) "plots/$filename"
                      else "$filename"
        }

        input:
        tuple val(samplename), path(bedgraph), path(reference_bacteria) from bedgraph_bact.groupTuple().combine(bact_table_len)

        output:
        path("*.html") into coverage_length_bacteria
        path("*_valid_bedgraph_files_bacteria") into valid_bedgraph_files_bacteria

        script:
        """
        generate_len_coverage_graph.py $samplename bacteria $reference_bacteria $bedgraph
        """
    }

    coverage_stats_bacteria.mix(failed_bacteria_samples).set { bacteria_coverage_results }

} else {

    nofile_path_bacteria_coverage = Channel.fromPath("NONE_bacteria")
    bacteria_results_template.combine(nofile_path_bacteria_coverage).set { bacteria_coverage_results }
}

if (params.fungi) {

    Channel.fromPath(params.fungi_dir_repo).into { fungi_table
                                                   fungi_table_len
                                                   fungi_sheet }

    if (params.fungi_ref_dir.endsWith('.gz') || params.fungi_ref_dir.endsWith('.tar') || params.fungi_ref_dir.endsWith('.tgz')) {

        process UNCOMPRESS_FUNGI_REF {
            label 'error_retry'

            input:
            path(ref_fungi) from params.fungi_ref_dir

            output:
            path("fungirefs") into fungi_ref_directory, fungi_references

            script:
            """
            mkdir "fungirefs"
            tar -xvf $ref_fungi --strip-components=1 -C "fungirefs"
            """
        }
    } else {
        Channel.fromPath("${params.fungi_ref_dir}").into { fungi_ref_directory
                                                           fungi_references }
    }

    process MASH_DETECT_FUNGI_REFERENCES {
        tag "$samplename"
        label "process_high"

        input:
        tuple val(samplename), val(single_end), path(reads), path(ref) from trimmed_fungi.combine(fungi_ref_directory)

        output:
        tuple val(samplename), path(mashout) into mash_result_fungi_references

        script:
        mashout = "mash_screen_results_fungi_${samplename}.txt"

        """
        find ${ref}/ -name "*" -type f > reference_list.txt
        mash sketch -k 32 -s 5000 -o reference -l reference_list.txt
        echo -e "#Identity\tShared_hashes\tMedian_multiplicity\tP-value\tQuery_id\tQuery_comment" > $mashout
        mash screen reference.msh $reads >> $mashout
        """
    }

    process SELECT_FINAL_FUNGI_REFERENCES {
        tag "$samplename"
        label "process_low"
        publishDir "${params.outdir}/${samplename}/fungi_coverage", mode: params.publish_dir_mode,
            saveAs { filename ->
                      if (filename.endsWith(".tsv")) filename
        }

        input:
        tuple val(samplename), path(mashresult), path(refdir), path(datasheet) from mash_result_fungi_references.combine(fungi_references).combine(fungi_sheet)

        output:
        tuple val(samplename), path("Final_fnas/*") optional true into bowtie_fungi_references
        tuple val(samplename), path("not_found.tsv") optional true into failed_fungi_samples

        script:
        """
        extract_significative_references.py $mashresult $refdir $datasheet

        """
    }

    trimmed_map_fungi.join(bowtie_fungi_references).set {bowtie_fungi_channel}

    def rawlist_fungi = bowtie_fungi_channel.toList().get()
    def bowtielist_fungi = []

    for (line in rawlist_fungi) {
        if (line[3] instanceof java.util.ArrayList){
            last_list = line[3]
            }
            else {
                last_list = [line[3]]
            }

            for (reference in last_list) {
                def ref_slice = [line[0],line[1],line[2],reference]
                bowtielist_fungi.add(ref_slice)
        }
    }

    def reads_fungi_mapping = Channel.fromList(bowtielist_fungi)

    process BOWTIE2_MAPPING_FUNGI {
        tag "${samplename} : ${reference}"
        label "process_high"

        input:
        tuple val(samplename), val(single_end), path(reads), path(reference) from reads_fungi_mapping

        output:
        tuple val(samplename), val(single_end), path("*_fungi.sam"), path(reference) into bowtie_alingment_sam_fungi

        script:
        samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"

        """
        bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $reference \\
        "index"

        bowtie2 \\
        -x "index" \\
        ${samplereads} \\
        -S "${reference}_vs_${samplename}_fungi.sam" \\
        --threads $task.cpus
        """
    }

    process SAMTOOLS_FUNGI {
        tag "$samplename"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(samfiles), path(reference) from bowtie_alingment_sam_fungi

        output:
        tuple val(samplename), val(single_end), path("*.sorted.bam") into bowtie_alingment_bam_fungi
        tuple val(samplename), val(single_end), path("*.sorted.bam"), path(reference) into ordered_bam_mpileup_fungi
        tuple val(samplename), val(single_end), path("*.sorted.bam.flagstat"), path("*.sorted.bam.idxstats"), path("*.sorted.bam.stats") into bam_stats_fungi

        script:

        """
        samtools view \\
        -@ $task.cpus \\
        -b \\
        -h \\
        -O BAM \\
        -o "\$(basename $samfiles .sam).bam" \\
        $samfiles

        samtools sort \\
        -@ $task.cpus \\
        -o "\$(basename $samfiles .sam).sorted.bam" \\
        "\$(basename $samfiles .sam).bam"

        samtools index "\$(basename $samfiles .sam).sorted.bam"

        samtools flagstat "\$(basename $samfiles .sam).sorted.bam" > "\$(basename $samfiles .sam).sorted.bam.flagstat"
        samtools idxstats "\$(basename $samfiles .sam).sorted.bam" > "\$(basename $samfiles .sam).sorted.bam.idxstats"
        samtools stats "\$(basename $samfiles .sam).sorted.bam" > "\$(basename $samfiles .sam).sorted.bam.stats"
        """
    }

    process BEDTOOLS_COVERAGE_FUNGI {
        tag "$samplename"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(bamfiles) from bowtie_alingment_bam_fungi

        output:
        tuple val(samplename), path("*_bedgraph.txt") into bedgraph_fungi
        tuple val(samplename), path("*_coverage.txt") into coverage_files_fungi_merge


        script:

        """
        bedtools genomecov -ibam $bamfiles  > "\$(basename -- $bamfiles .sorted.bam)_coverage.txt"
        bedtools genomecov -ibam $bamfiles -bga > "\$(basename -- $bamfiles .sorted.bam)_bedgraph.txt"
        """
    }

    process COVERAGE_STATS_FUNGI {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/fungi_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                     if (filename.endsWith(".html")) "plots/$filename"
                     else "$filename"

        }

        input:
        tuple val(samplename), path(coveragefiles), path(reference_fungi) from coverage_files_fungi_merge.groupTuple().combine(fungi_table)

        output:
        tuple val(samplename), path("*.tsv") into coverage_stats_fungi
        path("*.html") into coverage_graphs_fungi
        path("*_valid_coverage_files_fungi") into valid_coverage_files_fungi

        script:

        """
        graphs_coverage.py $samplename fungi $reference_fungi $coveragefiles
        """
    }

    process COVERAGE_LEN_FUNGI {

        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/fungi_coverage", mode: params.publish_dir_mode,
            saveAs: { filename ->
                      if (filename.endsWith(".html")) "plots/$filename"
                      }
        input:
        tuple val(samplename), path(bedgraph), path(reference_fungi) from bedgraph_bact.groupTuple().combine(fungi_table_len)

        output:
        path("*.html") into coverage_length_fungi
        path("*_valid_bedgraph_files_fungi") into valid_bedgraph_files_fungi

        script:
        """
        generate_len_coverage_graph.py $samplename fungi $reference_fungi $bedgraph
        """
    }

    coverage_stats_fungi.mix(failed_fungi_samples).set { fungi_coverage_results }

} else {

    nofile_path_fungi_coverage = Channel.fromPath("NONE_fungi")
    fungi_results_template.combine(nofile_path_fungi_coverage).set { fungi_coverage_results }
}

    if (params.translated_analysis) {

        process ASSEMBLY_METASPADES {
            tag "$samplename"
            label "process_high"
            publishDir "${params.outdir}/${samplename}/contigs", mode: params.publish_dir_mode

            input:
            tuple val(samplename), val(single_end), path(reads), val(filler) from reads_for_assembly

            output:
            tuple val(samplename), path("metaspades_result/contigs.fasta") into contigs, contigs_quast

            script:
            read = single_end ? "-s ${reads}" : "--meta -1 ${reads[0]} -2 ${reads[1]}"

            """
            spades.py \\
            $read \\
            --threads $task.cpus \\
            -o metaspades_result
            """
        }

        process QUAST_EVALUATION {
            tag "$samplename"
            label "process_medium"
            publishDir "${params.outdir}/${samplename}/quast_reports", mode: params.publish_dir_mode

            input:
            tuple val(samplename), file(contigfile) from contigs_quast

            output:
            file("$outputdir/report.html") into quast_results
            tuple val(samplename), path("$outputdir/report.tsv") into quast_multiqc
            path("$outputdir/report.tsv") into quast_multiqc_global

            script:
            outputdir = "quast_results_$samplename"

            """
            metaquast.py \\
            -f $contigfile \\
            -o $outputdir
            """
        }


        process KAIJU {
            tag "$samplename"
            label "process_high"

            input:
            tuple val(samplename), file(contig), path(kaijudb) from contigs.combine(kaiju_db)

            output:
            tuple val(samplename), path("*.out") into kaiju_results
            tuple val(samplename), path("*.krona") into kaiju_results_krona

            script:

            """
            kaiju \
            -t $kaijudb/nodes.dmp \\
            -f $kaijudb/*.fmi \\
            -i $contig \\
            -o ${samplename}_kaiju.out \\
            -z $task.cpus \\
            -v

            kaiju2table \\
            -t $kaijudb/nodes.dmp \\
            -n $kaijudb/names.dmp \\
            -r species \\
            -o ${samplename}_kaiju_summary.tsv \\
            ${samplename}_kaiju.out

            kaiju-addTaxonNames \\
            -t $kaijudb/nodes.dmp \\
            -n $kaijudb/names.dmp \\
            -i ${samplename}_kaiju.out \\
            -o ${samplename}_kaiju.names.out


            kaiju2krona \\
            -t $kaijudb/nodes.dmp \\
            -n $kaijudb/names.dmp \\
            -i ${samplename}_kaiju.out \\
            -o ${samplename}_kaiju.out.krona

            """
        }

        process KRONA_KAIJU_RESULTS {
            tag "$samplename"
            label "process_medium"
            publishDir "${params.outdir}/${samplename}/kaiju_results", mode: params.publish_dir_mode

            input:
            tuple val(samplename), path(kronafile) from kaiju_results_krona

            output:
            file("*.krona.html") into krona_results_kaiju

            script:
            outfile = "${samplename}_kaiju_result.krona.html"
            """
            ktImportText -o $outfile $kronafile
            """
        }

        process KAIJU_RESULTS_ANALYSIS {
            tag "$samplename"
            label "process_medium"
            publishDir "${params.outdir}/${samplename}/kaiju_results", mode: params.publish_dir_mode

            input:
            tuple val(samplename), path(outfile_kaiju) from kaiju_results

            output:
            tuple val(samplename), path("*_classified.txt"), path("*_unclassified.txt"), path("*_pieplot.html")

            script:
            """
            kaiju_results.py $samplename $outfile_kaiju
            """
        }

    }
    else {
        quast_multiqc_global = Channel.fromPath("NONE_quast_global")
        nofile_path_translated = Channel.fromPath("NONE_quast")
        samplechannel_translated.combine(nofile_path_translated).set { quast_multiqc }
    }

process GENERATE_RESULTS {
    tag "$samplename"
    label "process_low"
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(samplename), val(single_end), path(control_coverage), path(virus_coverage), path(bacteria_coverage), path(fungi_coverage) from sample_pe_template.join(control_coverage_results).join(virus_coverage_results).join(bacteria_coverage_results).join(fungi_coverage_results)

    output:
    path("*.html") into html_results
    val(samplename) into samplechannel_index

    script:
    paired = single_end ? "" : "--paired"
    control = (params.trimming && params.remove_control) ? "-control ${control_coverage}" : ""
    trimming = params.trimming ? "--trimming" : ""
    virus = params.virus ? "-virus '${virus_coverage}'" : ""
    bacteria = params.bacteria ? "-bacteria '${bacteria_coverage}'" : ""
    fungi = params.fungi ? "-fungi '${fungi_coverage}'" : ""
    scouting = params.kraken_scouting ? "--scouting" : ""
    translated_analysis = params.translated_analysis ? "--translated-analysis" : ""


    """
    generate-html.py \\
    --resultsdir ${params.outdir} \\
    --samplename $samplename \\
    $paired \\
    $control \\
    $trimming \\
    $virus \\
    $bacteria \\
    $fungi \\
    $scouting \\
    $translated_analysis

    """
}

/*
process MULTIQC_REPORT {
    tag "$samplename"
    label "process_medium"
    publishDir "${params.outdir}/${samplename}",  mode: params.publish_dir_mode

    input:
    tuple val(samplename), path(fastqc_raw), path(fastp_report), path(fastqc_trimmed), path(quast) from raw_fastqc_multiqc.join(fastp_multiqc).join(trimmed_fastqc_multiqc).join(quast_multiqc)

    output:
    path("*.html")

    script:

    """
    multiqc .
    """
}

process MULTIQC_REPORT_GLOBAL {

    label "process_medium"
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path(fastqc_raw) from raw_fastqc_multiqc_global.collect().ifEmpty([])
    path(fastp) from fastp_multiqc_global.collect().ifEmpty([])
    path(fastqc_trim) from trimmed_fastqc_multiqc_global.collect().ifEmpty([])
    path(quast) from quast_multiqc_global.collect().ifEmpty([])

    output:
    path("*.html")

    script:

    """
    multiqc .
    """
}
*/
process GENERATE_INDEX {
    label "process_low"
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    val(samplename_list) from samplechannel_index.collect()

    output:
    path("pikavirus_index.html") into pikavirus_index

    script:
    quality_control = params.trimming ? "--quality-control" : ""
    control_removal = params.remove_control ? "--control-removal" : ""
    scouting = params.kraken_scouting ? "--kraken_scouting" : ""

    virus = params.virus ? "--virus" : ""
    bacteria = params.bacteria ? "--bacteria" : ""
    fungi = params.fungi ? "--fungi" : ""

    translated_analysis = params.translated_analysis ? "--translated-analysis" : ""

    samplenames = samplename_list.join(" ")

    """
    create_index.py $quality_control \
                    $control_removal \
                    $scouting \
                    $virus \
                    $bacteria \
                    $fungi \
                    $translated_analysis \
                    --samplenames $samplenames
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[PikaVirus] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[PikaVirus] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.max_multiqc_email_size)
    // On success try attach the multiqc report

    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/pikavirus] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/pikavirus] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/pikavirus] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/pikavirus] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/pikavirus]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/pikavirus]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}

def isOffline() {
    try {
        return NXF_OFFLINE as Boolean
    }
    catch( Exception e ) {
        return false
    }
}

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
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
summary['Input']            = params.input
summary['Trimming']         = params.trimming
summary['Kaiju discovery']  = params.kaiju
if (params.kaiju) summary ['    Kaiju database']  = params.kaiju_db
summary['Virus Search']     = params.virus
if (params.virus) summary['    Virus Ref'] = params.vir_ref_dir
if (params.virus) summary['    Virus Index File'] = params.vir_dir_repo
summary['Bacteria Search']  = params.bacteria
if (params.bacteria) summary['    Bacteria Ref'] = params.bact_ref_dir
if (params.bacteria) summary['    Bacteria Index File'] = params.bact_dir_repo
summary['Fungi Search']     = params.fungi
if (params.fungi) summary['    Fungi Ref']     = params.fungi_ref_dir
if (params.fungi) summary['    Fungi Index File']     = params.fungi_dir_repo
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url

summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
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
    fastp -v &> v_fastp.txt &
    bowtie2 --version | head -n 1 | sed s/^.*version/version/ > v_bowtie2.txt
    mash -v | grep version &> v_mash.txt &
    spades.py -v &> v_spades.txt &
    quast -v &> v_quast.txt &
    samtools --version | grep samtools > v_samtools.txt
    bedtools -version > v_bedtools.txt
    
    kaiju -help &> tmp &
    head -n 1 tmp > v_kaiju.txt

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
if (params.kaiju){
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
    file "*_fastqc.{zip,html}" into fastqc_results
    tuple val(samplename), val(single_end), path("*.txt") into pre_filter_quality_data
    tuple val(samplename), path("*_fastqc.zip") into fastqc_multiqc_pre


    script:

    """
    fastqc --quiet --threads $task.cpus $reads

    for zipfile in *.zip;
    do
        unzip \$zipfile
        mv \$(basename \$zipfile .zip)/fastqc_data.txt \$(basename \$zipfile .zip).txt
    done
    """
}

/*
 * STEP 1.2 - TRIMMING​​​​​​​
*/
if (params.trimming) {
    process FASTP {
        tag "$samplename"
        label "process_medium"
        
        if (params.rescue_trimmed) {
            publishDir "${params.outdir}/${samplename}/trim_results", mode: params.publish_dir_mode,
        saveAs: { filename ->
                        filename.indexOf(".fastq") > 0 ? "trimmed/$filename" : "$filename"
                    }
        }

        input:
        tuple val(samplename), val(single_end), path(reads) from ch_cat_fortrim

        output:
        tuple val(samplename), val(single_end), path("*trim.fastq.gz") into trimmed_extract_virus, trimmed_extract_bact, trimmed_extract_fungi,
                                                                            reads_for_assembly, trimmed_paired_fastqc, trimmed_skip_hostremoval,
                                                                            trimmed_map_virus, trimmed_map_bact, trimmed_map_fungi
        
        tuple val(samplename), val(single_end), path("*fail.fastq.gz") into trimmed_unpaired

        script:
        detect_adapter =  single_end ? "" : "--detect_adapter_for_pe"
        reads1 = single_end ? "--in1 ${reads} --out1 ${samplename}_trim.fastq.gz --failed_out ${samplename}_fail.fastq.gz" : "--in1 ${reads[0]} --out1 ${samplename}_1_trim.fastq.gz --unpaired1 ${samplename}_1_fail.fastq.gz"
        reads2 = single_end ? "" : "--in2 ${reads[1]} --out2 ${samplename}_2_trim.fastq.gz --unpaired2 ${samplename}_2_fail.fastq.gz"
        
        """
        fastp \\
        $detect_adapter \\
        --cut_front \\
        --cut_tail \\
        --thread $task.cpus \\
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
        file "*_fastqc.{zip,html}" into trimmed_fastqc_results_html
        tuple val(samplename), path("*.txt") into post_filter_quality_data
        tuple val(samplename), path("*_fastqc.zip") into fastqc_multiqc_post

        script:
        
        """
        fastqc --quiet --threads $task.cpus $reads

        for zipfile in *.zip;
        do
            unzip \$zipfile
            mv \$(basename \$zipfile .zip)/fastqc_data.txt \$(basename \$zipfile .zip).txt
        done
        """
    }

    process EXTRACT_QUALITY_RESULTS {
        tag "$samplename"
        label "process_low"

        input:
        tuple val(samplename), val(single_end), path(pre_filter_data), path(post_filter_data) from pre_filter_quality_data.join(post_filter_quality_data)
        
        output:
        path("*_quality.txt") into quality_results_merged

        script:
        end = single_end ? "True" : "False"

        """
        extract_fastqc_data.py $samplename $params.outdir $end $pre_filter_data $post_filter_data

        """
    }


    process GENERATE_QUALITY_HTML {
        label "process_low"
        publishDir "${params.outdir}/quality_results", mode: params.publish_dir_mode

        input:
        path(quality_files) from quality_results_merged.collect()

        output:
        file("quality.html") into html_quality_result

        script:

        """
        cat $quality_files >> merged_file.txt
        
        merge_quality_stats.py merged_file.txt > quality.html

        """
    }
} else {
    Channel.from(ch_cat_fortrim).into{ trimmed_extract_virus trimmed_extract_bact trimmed_extract_fungi
                                      reads_for_assembly trimmed_paired_fastqc trimmed_skip_hostremoval
                                      trimmed_map_virus trimmed_map_bact trimmed_map_fungi }

}

if (params.virus) {    

    if (params.vir_ref_dir.endsWith('.gz') || params.vir_ref_dir.endsWith('.tar') || params.vir_ref_dir.endsWith('.tgz')) {

        process UNCOMPRESS_VIRUS_REF {
            label 'error_retry'

            input:
            path(ref_vir) from params.vir_ref_dir

            output:
            path("viralrefs") into virus_ref_directory, virus_references, virus_reference_graphcoverage
            
            script:
            """
            mkdir "viralrefs"
            tar -xvf $ref_vir --strip-components=1 -C "viralrefs"
            """
        }
    } else {
        Channel.fromPath("${params.vir_ref_dir}").into { virus_ref_directory
                                                         virus_references 
                                                         virus_reference_graphcoverage }
    }
   
    process MASH_DETECT_VIRUS_REFERENCES {
        tag "$samplename"
        label "process_high"
        
        input:
        tuple val(samplename), val(single_end), path(reads), path(ref) from trimmed_extract_virus.combine(virus_ref_directory)

        output:
        tuple val(samplename), path(mashout) into mash_result_virus_references

        script:
        mashout = "mash_results_virus_${samplename}.txt"
        
        """
        mash sketch -k 32 -s 5000 -r -m 2 -o query $reads
        mash sketch -k 32 -s 5000 -o reference ${ref}/*
        echo -e "#Reference-ID\tQuery-ID\tMash-distance\tP-value\tMatching-hashes" > $mashout

        mash dist reference.msh query.msh >> $mashout
        """       
    } 
    
    process SELECT_FINAL_VIRUS_REFERENCES {
        tag "$samplename"
        label "process_low"

        input:
        tuple val(samplename), path(mashresult), path(refdir) from mash_result_virus_references.combine(virus_references)

        output:
        tuple val(samplename), path("Final_fnas/*") into bowtie_virus_references

        script:
        """
        extract_significative_references.py $mashresult $refdir

        """
    }

    trimmed_map_virus.join(bowtie_virus_references).set{bowtie_virus_channel}

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

    def virus_reads_mapping = Channel.fromList(bowtielist_virus)

    process BOWTIE2_MAPPING_VIRUS {
        tag "${samplename} : ${reference}"
        label "process_high"
        
        input:
        tuple val(samplename), val(single_end), path(reads), path(reference) from virus_reads_mapping
        
        output:
        tuple val(samplename), val(single_end), path("*_virus.sam") into bowtie_alingment_sam_virus

        script:
        samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
        reference_name = reference.getFileName().toString()
        """
        bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $reference \\
        "index_${reference}"

        bowtie2 \\
        -x "index_${reference}" \\
        ${samplereads} \\
        -S "${reference}_vs_${samplename}_virus.sam" \\
        --threads $task.cpus
        
        """
    }

    process SAMTOOLS_BAM_FROM_SAM_VIRUS {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/virus_coverage/bam_stats", mode: params.publish_dir_mode

        input:
        tuple val(samplename), val(single_end), path(samfiles) from bowtie_alingment_sam_virus

        output:
        tuple val(samplename), val(single_end), path("*.sorted.bam") into bowtie_alingment_bam_virus
        tuple val(samplename), val(single_end), path("*.sorted.bam.flagstat"), path("*.sorted.bam.idxstats"), path("*.sorted.bam.stats") into bam_stats_virus
        
        script:

        """
        samtools view \\
        -@ $task.cpus \\
        -b \\
        -h \\
        -F4 \\
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

    process BEDTOOLS_COVERAGE_VIRUS {
        tag "$samplename"
        label "process_medium"

        input:
        tuple val(samplename), val(single_end), path(bamfiles) from bowtie_alingment_bam_virus

        output:
        tuple path("*_coverage_virus.txt"), path("*_bedgraph_virus.txt") into bedtools_coverage_files_virus
        tuple val(samplename), path("*_coverage_virus.txt") into coverage_files_virus_merge

        script:

        """
        bedtools genomecov -ibam $bamfiles -g "\$(basename -- $bamfiles .sorted.bam)_length.txt" > "\$(basename -- $bamfiles sorted.bam)_coverage_virus.txt"
        bedtools genomecov -ibam $bamfiles -g "\$(basename -- $bamfiles .sorted.bam)_length.txt" -bga >"\$(basename -- $bamfiles sorted.bam)_bedgraph_virus.txt"     
        """
    }
    
    process COVERAGE_STATS_VIRUS {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/virus_coverage", mode: params.publish_dir_mode

        input:
        tuple val(samplename), path(coveragefiles), path(reference_virus) from coverage_files_virus_merge.groupTuple().combine(virus_reference_graphcoverage)

        output:
        tuple val(samplename), path("*.csv") into coverage_stats_virus
        path("*.html") into coverage_graphs_virus
        
        script:
        outdirname = "${samplename}_virus"

        """
        graphs_coverage.py $outdirname $reference_virus $coveragefiles
        """        
    }
}

if (params.bacteria) {    

    if (params.bact_ref_dir.endsWith('.gz') || params.bact_ref_dir.endsWith('.tar') || params.bact_ref_dir.endsWith('.tgz')) {

        process UNCOMPRESS_BACTERIA_REF {
            label 'error_retry'

            input:
            path(ref_bact) from params.bact_ref_dir

            output:
            path("bactrefs") into bact_ref_directory, bact_references, bact_reference_graphcoverage
            
            script:
            """
            mkdir "bactrefs"
            tar -xvf $ref_bact --strip-components=1 -C "bactrefs"
            """
        }
    } else {
        Channel.fromPath("${params.bact_ref_dir}").into { bact_ref_directory bact_references bact_reference_graphcoverage  }
    }
   
    process MASH_DETECT_BACTERIA_REFERENCES {
        tag "$samplename"
        label "process_high"
        
        input:
        tuple val(samplename), val(single_end), path(reads), path(ref) from trimmed_extract_bact.combine(bact_ref_directory)

        output:
        tuple val(samplename), path(mashout) into mash_result_bact_references

        script:
        mashout = "mash_results_bact_${samplename}.txt"
        
        """
        mash sketch -k 32 -s 5000 -r -m 2 -o query $reads
        mash sketch -k 32 -s 5000 -o reference ${ref}/*
        echo -e "#Reference-ID\tQuery-ID\tMash-distance\tP-value\tMatching-hashes" > $mashout

        mash dist reference.msh query.msh >> $mashout
        """       
    } 
    
    process SELECT_FINAL_BACTERIA_REFERENCES {
        tag "$samplename"
        label "process_low"

        input:
        tuple val(samplename), path(mashresult), path(refdir) from mash_result_bact_references.combine(bact_references)

        output:
        tuple val(samplename), path("Final_fnas/*") into bowtie_bact_references

        script:
        """
        extract_significative_references.py $mashresult $refdir

        """
    }

    trimmed_map_bact.join(bowtie_bact_references).set{bowtie_bact_channel}

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

    def bact_reads_mapping = Channel.fromList(bowtielist_bact)

    process BOWTIE2_MAPPING_BACTERIA {
        tag "${samplename} : ${reference}"        
        label "process_high"
        
        input:
        tuple val(samplename), val(single_end), path(reads), path(reference) from bact_reads_mapping
        
        output:
        tuple val(samplename), val(single_end), path("*_bact.sam") into bowtie_alingment_sam_bact

        script:
        samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
        
        """
        bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $reference \\
        "index_${reference}"

        bowtie2 \\
        -x "index_${reference}" \\
        ${samplereads} \\
        -S "${reference}_vs_${samplename}_bact.sam" \\
        --threads $task.cpus
        """
    }

   process SAMTOOLS_BAM_FROM_SAM_BACTERIA {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/bacteria_coverage/bam_stats", mode: params.publish_dir_mode


        input:
        tuple val(samplename), val(single_end), path(samfiles) from bowtie_alingment_sam_bact

        output:
        tuple val(samplename), val(single_end), path("*.sorted.bam") into bowtie_alingment_bam_bact
        tuple val(samplename), val(single_end), path("*.sorted.bam.flagstat"), path("*.sorted.bam.idxstats"), path("*.sorted.bam.stats") into bam_stats_bact
        script:

        """
        samtools view \\
        -@ $task.cpus \\
        -b \\
        -h \\
        -F4 \\
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
        tuple path("*_coverage.txt"), path("*_bedgraph.txt") into bedtools_coverage_files_bact
        tuple val(samplename), path("*_coverage.txt") into coverage_files_bact_merge


        script:

        """
        bedtools genomecov -ibam $bamfiles -g "\$(basename -- $bamfiles .sorted.bam)_length.txt" > "\$(basename -- $bamfiles .sorted.bam)_coverage.txt"
        bedtools genomecov -ibam $bamfiles -g "\$(basename -- $bamfiles .sorted.bam)_length.txt" -bga >"\$(basename -- $bamfiles .sorted.bam)_bedgraph.txt"     
    
        """
    }
    
    process COVERAGE_STATS_BACTERIA {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/bacteria_coverage", mode: params.publish_dir_mode

        input:
        tuple val(samplename), path(coveragefiles), path(reference_bacteria) from coverage_files_bact_merge.groupTuple().combine(bact_reference_graphcoverage)

        output:
        tuple val(samplename), path("*.csv") into coverage_stats_bacteria
        path("*.html") into coverage_graphs_bacteria
        
        script:
        outdirname = "${samplename}_bacteria"

        """
        graphs_coverage.py $outdirname $reference_bacteria $coveragefiles
        """        
    }
}

if (params.fungi) {    

    if (params.fungi_ref_dir.endsWith('.gz') || params.fungi_ref_dir.endsWith('.tar') || params.fungi_ref_dir.endsWith('.tgz')) {

        process UNCOMPRESS_FUNGI_REF {
            label 'error_retry'

            input:
            path(ref_fungi) from params.fungi_ref_dir

            output:
            path("fungirefs") into fungi_ref_directory, fungi_references, fungi_reference_graphcoverage
            
            script:
            """
            mkdir "fungirefs"
            tar -xvf $ref_fungi --strip-components=1 -C "fungirefs"
            """
        }
    } else {
        Channel.fromPath("${params.fungi_ref_dir}").into { fungi_ref_directory fungi_references fungi_reference_graphcoverage }
    }
   
    process MASH_DETECT_FUNGI_REFERENCES {
        tag "$samplename"
        label "process_high"
        
        input:
        tuple val(samplename), val(single_end), path(reads), path(ref) from trimmed_extract_fungi.combine(fungi_ref_directory)

        output:
        tuple val(samplename), path(mashout) into mash_result_fungi_references

        script:
        mashout = "mash_results_fungi_${samplename}.txt"
        
        """
        mash sketch -k 32 -s 5000 -r -m 2 -o query $reads
        mash sketch -k 32 -s 5000 -o reference ${ref}/*
        echo -e "#Reference-ID\tQuery-ID\tMash-distance\tP-value\tMatching-hashes" > $mashout

        mash dist reference.msh query.msh >> $mashout
        """       
    } 
    
    process SELECT_FINAL_FUNGI_REFERENCES {
        tag "$samplename"
        label "process_low"

        input:
        tuple val(samplename), path(mashresult), path(refdir) from mash_result_fungi_references.combine(fungi_references)

        output:
        tuple val(samplename), path("Final_fnas/*") into bowtie_fungi_references

        script:
        """
        extract_significative_references.py $mashresult $refdir

        """
    }

    trimmed_map_fungi.join(bowtie_fungi_references).set{bowtie_fungi_channel}

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

    def fungi_reads_mapping = Channel.fromList(bowtielist_fungi)

    process BOWTIE2_MAPPING_FUNGI {
        tag "${samplename} : ${reference}"
        label "process_high"
        
        input:
        tuple val(samplename), val(single_end), path(reads), path(reference) from fungi_reads_mapping
        
        output:
        tuple val(samplename), val(single_end), path("*_fungi.sam") into bowtie_alingment_sam_fungi

        script:
        samplereads = single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
        
        """
        bowtie2-build \\
        --seed 1 \\
        --threads $task.cpus \\
        $reference \\
        "index_${reference}"

        bowtie2 \\
        -x "index_${reference}" \\
        ${samplereads} \\
        -S "${reference}_vs_${samplename}_fungi.sam" \\
        --threads $task.cpus
        """
    }

    process SAMTOOLS_BAM_FROM_SAM_FUNGI {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/fungi_coverage/bam_stats", mode: params.publish_dir_mode
        
        input:
        tuple val(samplename), val(single_end), path(samfiles) from bowtie_alingment_sam_fungi

        output:
        tuple val(samplename), val(single_end), path("*.sorted.bam") into bowtie_alingment_bam_fungi
        tuple val(samplename), val(single_end), path("*.sorted.bam.flagstat"), path("*.sorted.bam.idxstats"), path("*.sorted.bam.stats") into bam_stats_fungi
        script:

        """
        samtools view \\
        -@ $task.cpus \\
        -b \\
        -h \\
        -F4 \\
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
        tuple path("*_coverage.txt"), path("*_bedgraph.txt") into bedtools_coverage_files_fungi
        tuple val(samplename), path("*_coverage.txt") into coverage_files_fungi_merge


        script:

        """
        bedtools genomecov -ibam $bamfiles -g "\$(basename -- $bamfiles)_length.txt" > "\$(basename -- $bamfiles .sorted.bam)_coverage.txt"
        bedtools genomecov -ibam $bamfiles -g "\$(basename -- $bamfiles)_length.txt" -bga >"\$(basename -- $bamfiles .sorted.bam)_bedgraph.txt"        
        """
    }
    
    process COVERAGE_STATS_FUNGI {
        tag "$samplename"
        label "process_medium"
        publishDir "${params.outdir}/${samplename}/fungi_coverage", mode: params.publish_dir_mode

        input:
        tuple val(samplename), path(coveragefiles), path(reference_fungi) from coverage_files_fungi_merge.groupTuple().combine(fungi_reference_graphcoverage)

        output:
        tuple val(samplename), path("*.csv") into coverage_stats_fungi
        path("*.html") into coverage_graphs_fungi
        
        script:
        outdirname = "${samplename}_fungi"

        """
        graphs_coverage.py $outdirname $reference_fungi $coveragefiles
        """        
    }


}

if (params.translated_analysis) {

    if (params.kraken2_db.contains('.gz') || params.kraken2_db.contains('.tar')){

        process UNCOMPRESS_KRAKEN2DB {
            label 'error_retry'

            input:
            path(database) from params.kraken2_db

            output:
            path("kraken2db") into kraken2_db_files

            script:
            """
            mkdir "kraken2db"
            tar -zxf $database --strip-components=1 -C "kraken2db"
            """
        }
    } else {
        kraken2_db_files = Channel.fromPath(params.kraken2_db)
    }

    process SCOUT_KRAKEN2 {
        tag "$samplename"
        label "process_high"

        input: 
        tuple val(samplename), val(single_end), path(reads),path(kraken2db) from trimmed_paired_kraken2.combine(kraken2_db_files)

        output:
        tuple val(samplename), path("*.report") into kraken2_report_virus_references, kraken2_report_bacteria_references, kraken2_report_fungi_references
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
        trimmed_skip_hostremoval.into { reads_for_assembly }
    }

    if (params.kraken2krona) {

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
            file("*.krona.html") into krona_taxonomy

            script:
            outfile = "${samplename}_kraken.krona.html"

            """
            ktImportTaxonomy $kronafile -tax $taxonomy -o $outfile
            """
        }
    }

    if (params.kaiju){
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
            kaiju \\
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

}
/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/pikavirus] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/pikavirus] FAILED: $workflow.runName"
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
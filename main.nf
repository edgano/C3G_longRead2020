#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
params.bam = "$baseDir/data/RNA_nanopore.{brain,liver}.*.bam" 
//"http://www.genoscope.cns.fr/externe/ONT_mouse_RNA/data/transcriptome/RNA_nanopore.liver.C1R1_mapping_E94_minimap2_primary_no_read_less_than_80QC.bam"

// output directory
params.outdir = "${baseDir}/results"

log.info """\
         G  S  O  C    2   0    P  i  p  e  l  i  n  e  ~  version 0.1"
         =============================================================="
         Input sequences (BAM)                          : ${params.bam}
         Output directory (DIRECTORY)                   : ${params.outdir}
         """
         .stripIndent()


// Channels containing bam files
if ( params.bam ) {
  Channel
  .fromPath(params.bam)
  .map { item -> [ item.baseName.tokenize('.')[1], item] }
  .set { bam_ch }
}

process stringtie {
    tag "${pair_id}"
    publishDir "${params.outdir}/stringtie" , pattern: '*.gtf', mode: 'copy', overwrite: true
    container 'edgano/c3g'

    input:
      set val(pair_id), file(bam_file) from bam_ch

    output:
      set val(pair_id), file(bam_file), file("${pair_id}.out.gtf") into stringtie_ch

    script:
    """ 
    stringtie -L -o ${pair_id}.out.gtf ${bam_file}
    """
}
process tablemaker {
    tag "${pair_id}"
    publishDir "${params.outdir}/tablemaker",  pattern: '*.tablemaker', mode: 'copy', overwrite: true
    container 'edgano/c3g'

    input:
      set val(pair_id), file(bam_file), file(gtf_file) from stringtie_ch

    output:
      set val(pair_id), file(gtf_file),path("${pair_id}.tablemaker") into tablemaker_ch

    script:
    """
    tablemaker -W -G ${gtf_file} -o ${pair_id}.tablemaker ${bam_file}
    """
}

process ballgown {
    tag "${pair_id}"
    publishDir "${params.outdir}/ballgown" , pattern: '*_fpkm.csv', mode: 'copy', overwrite: true
    container 'kapeel/ballgown-r-package'

    input:
      set val(pair_id), file(gtf_file), file(tableFolder) from tablemaker_ch

    output:
      set val(pair_id), file(gtf_file), file("*_fpkm.csv") into ballgown_ch
    
    shell:
    '''
    #!/usr/bin/env Rscript

    library(ballgown)
    bg <- ballgown(dataDir = "!{baseDir}/results/tablemaker", samplePattern="!{pair_id}.tablemaker", meas='all')

    transcript_fpkm = texpr(bg, 'FPKM')
    write.csv(transcript_fpkm,file="!{pair_id}_fpkm.csv")
    '''
}

process postProduceOutputs{
    tag "${pair_id}"
    publishDir "${params.outdir}/gsoc" , mode: 'copy', overwrite: true

    input:
      set val(pair_id), file(gtf_file), file(fkm_file) from ballgown_ch

    output:
      set file("*.gtf"), file("*.csv") into result_ch

    script:
    """
    head -1000 ${gtf_file} > ${pair_id}_trim.gtf
    head -1000 ${fkm_file} > ${pair_id}_trim.csv
    """
}

workflow.onComplete {
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}"
}

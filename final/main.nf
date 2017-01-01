params.fqlist = "fqlist.txt"
params.fqpath = "$PWD/data/fq/"
params.genome = ""
params.gtf = "$PWD/data/genes.gtf"
params.fa = "$PWD/data/genome.fa"
params.p = 8
params.output = "results/"

/*
 * Initialization
 */

if( params.fqlist == "" ){
  error "--fqlist must be set."
}
else{
  fqlist_file = file(params.fqlist)
  if( !fqlist_file.exists() ){
    error "No fqlist file: ${fqlist_file}."
  }
  else{
    files_name = fqlist_file.text.split('\n|,')
    for( i = 0 ; i < files_name.length ; i = i + 1 ){
      current_file = file("${params.fqpath}${files_name[i]}")
      if( !current_file.exists() ){
        error "No fqfile file: ${current_file}."
      }
    }
  }
}

if( params.gtf == "" ){
  error "--gtf must be set."
}
else{
  gtf_file = file(params.gtf)
  if( !gtf_file.exists() ){
    error "No gtf file: ${gtf_file}."
  }
}

if( params.fa == "" ){
  error "--fa must be set."
}
else{
  fa_file = file(params.fa)
  if( !fa_file.exists() ){
    error "No fa file: ${fa_file}."
  }
}

if( params.genome == "" ){
  error "--genome must be set."
}
else{
  genome_folder = file(params.genome)
  if( !genome_folder.exists() ){
    error "No genome folder: ${genome_folder}."
  }
}

if( params.output == "" ){
  error "--output must be set."
}
else{
  output_folder = file(params.output)
  if( !output_folder.exists() ){
    output_folder.mkdirs()
  }
  assemblies_file = file("${output_folder}/assemblies.txt")
}

Channel
  .fromPath(fqlist_file)
  .splitCsv(header:["fq1","fq2"])
  .set { fq_pairs }

/*
 * Align the RNA-seq reads to the genome
 * 1| Map the reads for each sample to the reference genome.
 * Assemble expressed genes and transcripts
 * 2| Assemble transcripts for each sample.
 */

process AlignAndAssemble {

  maxForks 1

  input:
    set fq1, fq2 from fq_pairs

  output:
    set id, file('thout_folder'), file('clout_folder') into thcls
    file thout_folder
    file clout_folder

  exec:
    id = fq1.split('_')[0]
    fq1 = file("${params.fqpath}${fq1}")
    fq2 = file("${params.fqpath}${fq2}")

  script:
    """
    tophat -p $params.p -G $gtf_file -o thout_folder $genome_folder/genome $fq1 $fq2
    cufflinks -p $params.p -o clout_folder thout_folder/accepted_hits.bam
    """
}

/*
 * 3| Create a file called assemblies.txt that lists the assembly file for each sample.
 */
process CreateAssemblyFile{
  input:
    set id, thout_folder, clout_folder from thcls

  output:
    file assembly_file

  exec:
    thout_folder.moveTo("${output_folder}/${id}_thout")
    clout_folder.moveTo("${output_folder}/${id}_clout")

  script:
    """
    touch assemblies_file
    echo ${output_folder}/${id}_clout/transcripts.gtf > assembly_file
    """
}

assemblies_file = assembly_file.collectFile( name:"${output_folder}/assemblies.txt")

/*
 * 4| Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation.
 */

process Cuffmerge {

  input:
    file assemblies_file

  output:
    stdout channel

  script:
    """
    cuffmerge -g $gtf_file -s $fa_file -p 8 ${assemblies_file}
    """

}

channel.subscribe{
  println it
}

/*
 * Identify differentially expressed genes and transcripts
 * 5| Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate.
 */
//
// process CuffdiffAnd {
//
//   input:
//
//   output:
//
//   exec:
//
//   script:
//
// }
//
// Align the RNA-seq reads to the genome
// 1| Map the reads for each sample to the reference genome.
// $ tophat -p 8 -G genes.gtf -o C1_R1_thout genome C1_R1_1.fq C1_R1_2.fq
// $ tophat -p 8 -G genes.gtf -o C1_R2_thout genome C1_R2_1.fq C1_R2_2.fq
// $ tophat -p 8 -G genes.gtf -o C1_R3_thout genome C1_R3_1.fq C1_R3_2.fq
// $ tophat -p 8 -G genes.gtf -o C2_R1_thout genome C2_R1_1.fq C1_R1_2.fq
// $ tophat -p 8 -G genes.gtf -o C2_R2_thout genome C2_R2_1.fq C1_R2_2.fq
// $ tophat -p 8 -G genes.gtf -o C2_R3_thout genome C2_R3_1.fq C1_R3_2.fq
//
// Assemble expressed genes and transcripts
// 2| Assemble transcripts for each sample:
// $ cufflinks -p 8 -o C1_R1_clout C1_R1_thout/accepted_hits.bam
// $ cufflinks -p 8 -o C1_R2_clout C1_R2_thout/accepted_hits.bam
// $ cufflinks -p 8 -o C1_R3_clout C1_R3_thout/accepted_hits.bam
// $ cufflinks -p 8 -o C2_R1_clout C2_R1_thout/accepted_hits.bam
// $ cufflinks -p 8 -o C2_R2_clout C2_R2_thout/accepted_hits.bam
// $ cufflinks -p 8 -o C2_R3_clout C2_R3_thout/accepted_hits.bam
//
// 3| Create a file called assemblies.txt that lists the assembly file for each sample. The file should contain the following lines:
// ./C1_R1_clout/transcripts.gtf
// ./C2_R2_clout/transcripts.gtf
// ./C1_R2_clout/transcripts.gtf
// ./C2_R1_clout/transcripts.gtf
// ./C1_R3_clout/transcripts.gtf
// ./C2_R3_clout/transcripts.gtf
//
// 4| Run Cuffmerge on all your assemblies to create a single merged transcriptome annotation:
// cuffmerge -g genes.gtf -s genome.fa -p 8 assemblies.txt
//
// Identify differentially expressed genes and transcripts ● TIMING ~6 h
// 5| Run Cuffdiff by using the merged transcriptome assembly along with the BAM files from TopHat for each replicate:
// $ cuffdiff -o diff_out -b genome.fa -p 8 –L C1,C2 -u merged_asm/merged.gtf \
// ./C1_R1_thout/accepted_hits.bam,./C1_R2_thout/accepted_hits.bam,./C1_R3_thout/
// accepted_hits.bam \
// ./C2_R1_thout/accepted_hits.bam,./C2_R3_thout/accepted_hits.bam,./C2_R2_thout/
// accepted_hits.bam

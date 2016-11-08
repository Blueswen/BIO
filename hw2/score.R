library('proto')
library('argparse')
library('seqinr')
Rstudio <- FALSE

#####################
## Process options ##
#####################

# Rscript score.R --input result.fasta --score PAM250.txt --gap_open -10 --gap_extend -2

parser <- ArgumentParser(description='Calculate pairwise sequence alignment score.')
parser$add_argument('--input', help='input fasta file', required=TRUE)
parser$add_argument('--score', help='input score file', required=TRUE)
parser$add_argument('--gap_open', type='integer', help='gap open penalty', default=-10)
parser$add_argument('--gap_extend', type='integer', help='gap extend penalty', default=-2)
if(Rstudio){
  m_input <- '~/Documents/NCCU/1051/Bio/hw/hw2/result.fasta'
  m_score <- '~/Documents/NCCU/1051/Bio/hw/hw2/PAM250.txt'
  m_gap_open <- -10
  m_gap_extend <- -2
  args <- parser$parse_args(c('--input',m_input,'--score',m_score,'--gap_open',m_gap_open,'--gap_extend',m_gap_extend))
}else{
  args <- parser$parse_args()
}

input <- args$input
score <- args$score
gap_open <- args$gap_open
gap_extend <- args$gap_extend

#####################
##    Load Data    ##
#####################

input_fasta <-read.alignment(file=input, format='fasta', forceToLower=FALSE)
if(input_fasta$nb<2){
  stop("input_file must contains more than 2 sequences.")
}else if(input_fasta$nb>2){
  print("There are more than 2 sequneces in input_file, first and second sequence will be the targets.")
}else if(nchar(input_fasta$seq[1])!=nchar(input_fasta$seq[2])){
  stop("First 2 sequences length are not equal.")
}

score_table <- read.table(file=score)
if( nrow(score_table)!=24 | ncol(score_table)!=24 ){
  stop("score_file format isn't correct.")
}

#####################
##       Main      ##
#####################

seq1 <- input_fasta$seq[1]
seq2 <- input_fasta$seq[2]

open_flag = FALSE
align_score = 0
for( i in 1:nchar(seq1)){
  ac1 = substr(seq1,i,i)
  ac2 = substr(seq2,i,i)
  if( ac1!='-' & ac2!='-' ){
    if(open_flag){
      open_flag = FALSE
    }
    align_score = align_score + score_table[ac1,ac2]
  }else{
    if(open_flag){
      align_score = align_score + gap_extend
    }else{
      open_flag = TRUE
      align_score = align_score + gap_open
    }
  }
}

print(paste('Score: ',align_score))

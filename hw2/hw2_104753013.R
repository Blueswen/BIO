library('proto')
library('argparse')
library('seqinr')
Rstudio <- FALSE

#####################
## Process options ##
#####################

# Rscript hw2_104753013.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta

args <- NULL
parser <- ArgumentParser(description='Pairwise sequence alignment with gap open penalty and gap extend penalty.')
parser$add_argument('--input', help='input fasta file', required=TRUE) 
parser$add_argument('--score', help='input score file', required=TRUE)
parser$add_argument('--aln', type='character', choices=c('global','local'), default='global', help='align way')
parser$add_argument('--gap_open', type='double', help='gap open num', default=-10)
parser$add_argument('--gap_extend', type='double', help='gap extend num', default=-2)
parser$add_argument('--output', help='output fasta file', required=TRUE)
parser$add_argument('--no_penalty', action='store_true', help='disable two gap-penalty scheme')

if(Rstudio){
  m_input <- '~/Documents/NCCU/1051/Bio/hw/hw2/test.fasta'
  m_score <- '~/Documents/NCCU/1051/Bio/hw/hw2/pam250.txt'
  m_aln <- 'local'
  m_gap_open <- -10
  m_gap_extend <- -2
  m_output <- '~/Documents/NCCU/1051/Bio/hw/hw2/result.fasta'
  args <- parser$parse_args(c('--input',m_input,'--score',m_score,'--aln',m_aln,'--gap_open',m_gap_open,'--gap_extend',m_gap_extend,'--output',m_output,'--no_penalty'))
}else{
  args <- parser$parse_args()
}
   
input <- args$input
score <- args$score
aln <- args$aln
gap_open <- args$gap_open
gap_extend <- args$gap_extend
output <- args$output
penalty <- args$no_penalty

if(!file.exists(input)){
  stop("input_file doesn't exist.")
}

if(!file.exists(score)){
  stop("score_file doesn't exist.")
}

if(file.exists(output)){
  stop("output_file already exists.")
}


#####################
##    Load Data    ##
#####################
input_fasta <-read.alignment(file=input, format='fasta', forceToLower=FALSE)
if(input_fasta$nb<2){
  stop("input_file must contains more than 2 sequences.")
}else if(input_fasta$nb>2){
  print("There are more than 2 sequneces in input_file, first and second sequence will be the targets.")
}

score_table <- read.table(file=score)
if( nrow(score_table)!=24 | ncol(score_table)!=24 ){
  print("score_file format isn't correct, it could cause error.")
}

#####################
##       Main      ##
#####################

fa <- c(input_fasta$seq[1], input_fasta$seq[2])

# Initial table
main_table <- data.frame(matrix(ncol = nchar(fa[1])+1, nrow = nchar(fa[2])+1))
main_col = append('-',unlist(strsplit(fa[1],"{1}")))
main_row = append('-',unlist(strsplit(fa[2],"{1}")))

Ix_table <- data.frame(matrix(ncol = nchar(fa[1])+1, nrow = nchar(fa[2])+1))
Iy_table <- data.frame(matrix(ncol = nchar(fa[1])+1, nrow = nchar(fa[2])+1))

alignments = c("","")
cur_i = 0
cur_j = 0

# Fill main table
if(aln == "global" & penalty == FALSE){
  # Assign first row and first col
  for( i in 1:length(main_row)){
    #main_table[i,1] = (i-1)*gap_open
    main_table[i,1] = 0
  }
  for( j in 1:length(main_col)){
    #main_table[1,j] = (j-1)*gap_open
    main_table[1,j] = 0
  }
  
  for( i in 2:length(main_row)){
    for( j in 2:length(main_col)){
      sub_s = score_table[main_row[i], main_col[j]] + main_table[i-1,j-1]
      del_s = gap_open + main_table[i-1,j]
      ins_s = gap_open + main_table[i,j-1]
      main_table[i,j] = max(c(sub_s, del_s, ins_s))
    }  
  }
  
  cur_i = length(main_row)
  cur_j = length(main_col)
  
  # Trace-back
  while( cur_i>1 & cur_j>1 ){
    if(main_table[cur_i,cur_j] == main_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i],main_col[cur_j]]){
      alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
      alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
      cur_i = cur_i - 1
      cur_j = cur_j - 1
    }else if(main_table[cur_i,cur_j] == main_table[cur_i-1,cur_j] + gap_open){
      alignments[1] = paste('-', alignments[1], sep="")
      alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
      cur_i = cur_i - 1
    }else if(main_table[cur_i,cur_j] == main_table[cur_i,cur_j-1] + gap_open){
      alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
      alignments[2] = paste('-', alignments[2], sep="")
      cur_j = cur_j - 1
    }
  }
  while( cur_i>1 ){
    alignments[1] = paste('-', alignments[1], sep="")
    alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
    cur_i = cur_i - 1
  }
  while( cur_j>1 ){
    alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
    alignments[2] = paste('-', alignments[2], sep="")
    cur_j = cur_j - 1
  }
  
}else if(aln == "local" & penalty == FALSE){
  # Assign first row and first col
  for( i in 1:length(main_row)){
    main_table[i,1] = 0
  }
  for( j in 1:length(main_col)){
    main_table[1,j] = 0
  }
  
  best_score = 0
  best_i = 0
  best_j = 0
  
  for( i in 2:length(main_row)){
    for( j in 2:length(main_col)){
      sub_s = score_table[main_row[i], main_col[j]] + main_table[i-1,j-1]
      del_s = gap_open + main_table[i-1,j]
      ins_s = gap_open + main_table[i,j-1]
      main_table[i,j] = max(c(sub_s, del_s, ins_s, 0))
      
      #prepare for trace back
      if(main_table[i,j] >= best_score){
        best_score = main_table[i,j]
        best_i = i
        best_j = j
      }
    }  
  }

  cur_i = best_i
  cur_j = best_j
  
  # Trace-back
  while( cur_i>1 & cur_j>1 ){
    if(main_table[cur_i,cur_j] == main_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i],main_col[cur_j]]){
      alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
      alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
      cur_i = cur_i - 1
      cur_j = cur_j - 1
    }else if(main_table[cur_i,cur_j] == main_table[cur_i-1,cur_j] + gap_open){
      alignments[1] = paste('-', alignments[1], sep="")
      alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
      cur_i = cur_i - 1
    }else if(main_table[cur_i,cur_j] == main_table[cur_i,cur_j-1] + gap_open){
      alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
      alignments[2] = paste('-', alignments[2], sep="")
      cur_j = cur_j - 1
    }
  }
  
}else if(aln == 'global' & penalty == TRUE){
  # Assign first row and first col
  main_table[1,1] = 0
  Ix_table[1,1] = 0
  Iy_table[1,1] = 0
  for( i in 2:length(main_row)){
    main_table[i,1] = gap_open + (i-2)*gap_extend
    Ix_table[i,1] = 0
    Iy_table[i,1] = -Inf
  }
  for( j in 2:length(main_col)){
    main_table[1,j] = gap_open + (j-2)*gap_extend
    Ix_table[1,j] = -Inf
    Iy_table[1,j] = 0
  }
  
  for( i in 2:length(main_row)){
    for( j in 2:length(main_col)){
      # Ix_table
      Ix_open_s = main_table[i-1,j] + gap_open
      Ix_extend_s = Ix_table[i-1,j] + gap_extend
      Ix_table[i,j] = max(c(Ix_open_s,Ix_extend_s))
      # Iy_table
      Iy_open_s = main_table[i,j-1] + gap_open
      Iy_extend_s = Iy_table[i,j-1] + gap_extend
      Iy_table[i,j] = max(c(Iy_open_s,Iy_extend_s))
      # main_table
      sub_s = main_table[i-1,j-1] + score_table[main_row[i], main_col[j]]
      Ix_s = Ix_table[i-1,j-1] + score_table[main_row[i], main_col[j]]
      Iy_s = Iy_table[i-1,j-1] + score_table[main_row[i], main_col[j]]
      main_table[i,j] = max(c(sub_s, Ix_s, Iy_s))
    }  
  }
  cur_i = length(main_row)
  cur_j = length(main_col)
  cur_score = max(main_table[cur_i,cur_j], Ix_table[cur_i,cur_j], Iy_table[cur_i,cur_j])
  table_ind = which.max(c(main_table[cur_i,cur_j], Ix_table[cur_i,cur_j], Iy_table[cur_i,cur_j]))
  
  # Trace-back
  while( cur_i>1 & cur_j>1 ){
    if(table_ind == 1){
      if(cur_score == main_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i],main_col[cur_j]]){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_j = cur_j - 1
        cur_score = main_table[cur_i,cur_j]
        table_ind = 1
      }else if(cur_score == Ix_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i], main_col[cur_j]]){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_j = cur_j - 1
        cur_score = Ix_table[cur_i,cur_j]
        table_ind = 2
      }else if(cur_score == Iy_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i], main_col[cur_j]]){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_j = cur_j - 1
        cur_score = Iy_table[cur_i,cur_j]
        table_ind = 3
      }
    }else if(table_ind == 2){
      if(cur_score == main_table[cur_i-1,cur_j] + gap_open){
        alignments[1] = paste('-', alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_score = main_table[cur_i,cur_j]
        table_ind = 1
      }else if(cur_score == Ix_table[cur_i-1,cur_j] + gap_extend){
        alignments[1] = paste('-', alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_score = Ix_table[cur_i,cur_j]
      }
    }else if(table_ind == 3){
      if(cur_score == main_table[cur_i,cur_j-1] + gap_open){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste('-', alignments[2], sep="")
        cur_j = cur_j - 1
        cur_score = main_table[cur_i,cur_j]
        table_ind = 1
      }else if(cur_score == Iy_table[cur_i,cur_j-1] + gap_extend){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste('-', alignments[2], sep="")
        cur_j = cur_j - 1
        cur_score = Iy_table[cur_i,cur_j]
      }
    }
  }
  while( cur_i>1 ){
    alignments[1] = paste('-', alignments[1], sep="")
    alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
    cur_i = cur_i - 1
  }
  while( cur_j>1 ){
    alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
    alignments[2] = paste('-', alignments[2], sep="")
    cur_j = cur_j - 1
  }
  
}else if(aln == 'local' & penalty == TRUE){
  best_score = 0
  best_i = 0
  best_j = 0
  
  # Assign first row and first col
  main_table[1,1] = 0
  Ix_table[1,1] = 0
  Iy_table[1,1] = 0
  for( i in 2:length(main_row)){
    main_table[i,1] = gap_open + (i-2)*gap_extend
    Ix_table[i,1] = 0
    Iy_table[i,1] = -Inf
  }
  for( j in 2:length(main_col)){
    main_table[1,j] = gap_open + (j-2)*gap_extend
    Ix_table[1,j] = -Inf
    Iy_table[1,j] = 0
  }
  
  for( i in 2:length(main_row)){
    for( j in 2:length(main_col)){
      # Ix_table
      Ix_open_s = main_table[i-1,j] + gap_open
      Ix_extend_s = Ix_table[i-1,j] + gap_extend
      Ix_table[i,j] = max(c(Ix_open_s,Ix_extend_s,0))
      # Iy_table
      Iy_open_s = main_table[i,j-1] + gap_open
      Iy_extend_s = Iy_table[i,j-1] + gap_extend
      Iy_table[i,j] = max(c(Iy_open_s,Iy_extend_s,0))
      # main_table
      sub_s = main_table[i-1,j-1] + score_table[main_row[i], main_col[j]]
      Ix_s = Ix_table[i-1,j-1] + score_table[main_row[i], main_col[j]]
      Iy_s = Iy_table[i-1,j-1] + score_table[main_row[i], main_col[j]]
      main_table[i,j] = max(c(sub_s, Ix_s, Iy_s,0))
      
      #prepare for trace back
      if(max(c(main_table[i,j], Ix_table[i,j], Iy_table[i,j])) >= best_score){
        best_score = max(c(main_table[i,j], Ix_table[i,j], Iy_table[i,j]))
        best_i = i
        best_j = j
      }
    } 
  }
  
  cur_i = best_i
  cur_j = best_j
  cur_score = max(c(main_table[cur_i,cur_j], Ix_table[cur_i,cur_j], Iy_table[cur_i,cur_j]))
  table_ind = which.max(c(main_table[cur_i,cur_j], Ix_table[cur_i,cur_j], Iy_table[cur_i,cur_j]))
  
  # Trace-back
  while( cur_i>1 & cur_j>1 ){
    if(table_ind == 1){
      if(cur_score == main_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i],main_col[cur_j]]){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_j = cur_j - 1
        cur_score = main_table[cur_i,cur_j]
        table_ind = 1
      }else if(cur_score == Ix_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i], main_col[cur_j]]){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_j = cur_j - 1
        cur_score = Ix_table[cur_i,cur_j]
        table_ind = 2
      }else if(cur_score == Iy_table[cur_i-1,cur_j-1] + score_table[main_row[cur_i], main_col[cur_j]]){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_j = cur_j - 1
        cur_score = Iy_table[cur_i,cur_j]
        table_ind = 3
      }
    }else if(table_ind == 2){
      if(cur_score == main_table[cur_i-1,cur_j] + gap_open){
        alignments[1] = paste('-', alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_score = main_table[cur_i,cur_j]
        table_ind = 1
      }else if(cur_score == Ix_table[cur_i-1,cur_j] + gap_extend){
        alignments[1] = paste('-', alignments[1], sep="")
        alignments[2] = paste(main_row[cur_i], alignments[2], sep="")
        cur_i = cur_i - 1
        cur_score = Ix_table[cur_i,cur_j]
      }
    }else if(table_ind == 3){
      if(cur_score == main_table[cur_i,cur_j-1] + gap_open){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste('-', alignments[2], sep="")
        cur_j = cur_j - 1
        cur_score = main_table[cur_i,cur_j]
        table_ind = 1
      }else if(cur_score == Iy_table[cur_i,cur_j-1] + gap_extend){
        alignments[1] = paste(main_col[cur_j], alignments[1], sep="")
        alignments[2] = paste('-', alignments[2], sep="")
        cur_j = cur_j - 1
        cur_score = Iy_table[cur_i,cur_j]
      }
    }
  }
  
}

#####################
##      Output     ##
#####################

if(!Rstudio){
  write.fasta(as.list(alignments),input_fasta$nam,output)
}

# Calculate Score
open_flag = FALSE
align_score = 0
for( i in 1:nchar(alignments[1])){
  ac1 = substr(alignments[1],i,i)
  ac2 = substr(alignments[2],i,i)
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
print('Alignment result:')
print(alignments[1])
print(alignments[2])
print(paste('Score: ',align_score))
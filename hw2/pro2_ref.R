######################################
# the reference code of program2 
######################################

######################################
# initial
######################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<-args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-args[i+1]
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

print("PARAMETERS")
print(paste("input file         :", i_f))
print(paste("output file        :", o_f))
print(paste("score file         :", s_f))
print(paste("gap open penalty   :", g_o))
print(paste("gap extend penalty :", g_e))

######################################
# main
######################################
# your code
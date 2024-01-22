source("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/node_summary_funs.R")

args <- commandArgs(T)

rds <- args[1]
social <- args[2]
subsocial <- args[3]

seq_out <- node_seq_name(rds,social,subsocial)
write_tsv(seq_out,paste(social,"_",subsocial,"_uniq_node_seq.tsv",sep=""))

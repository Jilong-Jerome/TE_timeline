from gwf import *
from workflow_templates import *
from workflow_targets import *


gwf = Workflow()

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/Ste_ncbi_chromosome.fa"
name = "stegodyphus"
make_db(gwf,path,fasta,name)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum"
seed = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/N1168_seq.fa"
outname = "N1168"
db = "stegodyphus"
run_blast_search(gwf,path,db,seed,outname)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum"
seed = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum/N1168_blast_pass_dedup.fa"
outname = "N1168_blast_r2"
db = "stegodyphus"
run_blast_search(gwf,path,db,seed,outname)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum"
seed = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum/N1168_blast_sar.fa"
outname = "N1168_sar_r2"
db = "stegodyphus"
run_blast_search(gwf,path,db,seed,outname)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum"
seed = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/N305_seq.fa"
outname = "N305"
db = "stegodyphus"
run_blast_search(gwf,path,db,seed,outname)
## Prepare pair uniq nodes
social = "dum"
subsocial = "tent"
rds = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/graphs/dum_tent_wgs_comp2_self_blast_sim95.rds"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/dum_tent"
run_find_uniq_seq(gwf,path,rds,social,subsocial)
run_find_uniq_seq2fa(gwf,path,social,subsocial)
run_seq_per_node(gwf,path,social,subsocial)

social = "sar"
subsocial = "bic"
rds = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/graphs/sar_bic_wgs_comp2_self_blast_sim95.rds"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/sar_bic"
run_find_uniq_seq(gwf,path,rds,social,subsocial)
run_find_uniq_seq2fa(gwf,path,social,subsocial)
run_seq_per_node(gwf,path,social,subsocial)

# blast run per node
pair = "dum_tent"
for sp in ["dum","tent"]:
    nodes = open("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/dum_tent/{sp}_node.txt".format(sp=sp))
    for node in nodes:
        node = node.strip("\n")
        run_build_node_seed(gwf,pair,sp,node)
        if sp == "tent":
            new_sp = "ten"
        else:
            new_sp = sp
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}".format(sp=new_sp,node=node)
        seed = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{node}_seq.fa".format(sp=new_sp,node=node)
        outname = new_sp+"_"+node
        db = "stegodyphus"
        run_blast_search_node(gwf,path,db,seed,node,outname)
        run_blast_compose(gwf,new_sp,node)
        run_blast_align(gwf,new_sp,node)
        run_blast_dist(gwf,new_sp,node)
        run_node_root(gwf,new_sp,node)
        run_node_timeline(gwf,new_sp,node)

# Node analysis
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum/N1168_raxml"
aln = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/dum_tent/dum/N1168_blast_pass_dedup_aln_filtered.fasta"
outname = "N1168_ML"
run_raxml_build(gwf,path,aln,outname)

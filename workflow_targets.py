from workflow_templates import *

def make_db(gwf,path,fasta,name):
    gwf.target_from_template(
    name = "blastdb_{name}".format(name = name),
    template = blast_makedb(path=path,fasta=fasta,name=name)
    )

def run_blast_search(gwf,path,db,seed,outname):
    gwf.target_from_template(
    name = "blast_{outname}".format(outname=outname),
    template = blast_search(path,db,seed,outname)
    )

def run_blast_search_node(gwf,path,db,seed,node,outname):
    gwf.target_from_template(
    name = "blast_{outname}".format(outname=outname),
    template = blast_search_node(path,db,seed,node,outname)
    )

def run_find_uniq_seq(gwf,path,rds,social,subsocial):
    gwf.target_from_template(
    name = "rds2node_summary_{social}_{subsocial}".format(social=social,subsocial=subsocial),
    template = graph2uniq_node_seq(path,rds,social,subsocial)
    )
def run_find_uniq_seq2fa(gwf,path,social,subsocial):
    gwf.target_from_template(
    name = "uniqnode2fa_{social}_{subsocial}".format(social=social,subsocial=subsocial),
    template = unique_node_seq(path,social,subsocial)
    )
def run_seq_per_node(gwf,path,social,subsocial):
    seq_summary = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/{social}_{subsocial}/{social}_{subsocial}_uniq_node_seq.tsv".format(social=social,subsocial=subsocial)
    gwf.target_from_template(
    name = "seq_per_node_{social}_{subsocial}".format(social=social,subsocial=subsocial),
    template = make_node_txt(path,social,subsocial,seq_summary)
    )

def run_build_node_seed(gwf,pair,sp,node):
    if sp == "tent":
        new_sp = "ten"
    else:
        new_sp = sp
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}".format(sp=new_sp,node=node)
    gwf.target_from_template(
    name = "build_blast_seed_{sp}_{node}".format(sp=sp,node=node),
    template = build_node_seed(path,pair,sp,node)
    )

def run_blast_compose(gwf,sp,node):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}".format(sp=sp,node=node)
    gwf.target_from_template(
    name = "blast_compose_{sp}_{node}".format(sp=sp,node=node),
    template = blast_compose(path,sp,node)
    )

def run_blast_align(gwf,sp,node):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}".format(sp=sp,node=node)
    gwf.target_from_template(
    name = "blast_align_{sp}_{node}".format(sp=sp,node=node),
    template = align_blast(path,sp,node)
    )

def run_blast_dist(gwf,sp,node):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_2_phylo/{sp}/{node}".format(sp=sp,node=node)
    gwf.target_from_template(
    name = "blast_dist_{sp}_{node}".format(sp=sp,node=node),
    template = build_aln_distance(path,sp,node)
    )

def run_node_root(gwf,sp,node):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_2_phylo/{sp}/{node}".format(sp=sp,node=node)
    gwf.target_from_template(
    name = "set_root_{sp}_{node}".format(sp=sp,node=node),
    template = root_node(path,sp,node)
    )

def run_node_timeline(gwf,sp,node):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_3_timeline/{sp}/{node}/data".format(sp=sp,node=node)
    if sp == "dum":
        sister_sp = "ten"
    elif sp == "ten" or sp == "tent":
        sister_sp = "dum"
    elif sp == "sar":
        sister_sp = "bic"
    elif sp == "bic":
        sister_sp = "sar"
    gwf.target_from_template(
    name = "timeline_{sp}_{node}".format(sp=sp,node=node),
    template = node_timeline(path,sp,sister_sp,node)
    )

def run_raxml_build(gwf,path,aln,outname):
    gwf.target_from_template(
    name = "raxml_{outname}".format(outname = outname),
    template = raxml_tree(path,aln,outname) 
    )

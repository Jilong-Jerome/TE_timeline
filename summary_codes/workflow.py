from gwf import *
from workflow_templates import *
from workflow_targets import *

SCRIPT = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/summary_codes"
ROOT = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_4_summary"
LOG = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/logs/orf"

gwf = Workflow()


# Prepare
quantile = 3 # quantile% for defining speciation time as the scaler of branch length across different nodes 
for pair in ["dum_ten"]:
    ind_nodes = open(SCRIPT+ "/{pair}_independent_node_list.tsv".format(pair=pair))
    for node in ind_nodes:
        infos = node.strip("\n").split("\t")
        node = infos[0]
        sp = infos[1]
        condition = infos[2]
        if node != "node":
            prepare_nodes_to_summary(gwf,ROOT,LOG,sp,node,pair,condition)
            for quantile in [0.5,1,3,5]:
                cal_scaled_timeline(gwf,ROOT,LOG,SCRIPT,sp,node,pair,condition,quantile)
                cal_scaled_terminal(gwf,ROOT,LOG,SCRIPT,sp,node,pair,condition,quantile)

from workflow_templates import *

def prepare_nodes_to_summary(gwf,root_path,log_path,sp,node,pair,condition):
    path = root_path +  "/{pair}/{condition}/{node}".format(pair=pair,condition=condition,node=node)
    gwf.target_from_template(
    name = "prepare_independent_{node}".format(node=node),
    template = retrieve_independent_nodes(path,log_path,sp,node)
    )

def cal_scaled_timeline(gwf,root_path,log_path,script_path,sp,node,pair,condition,quantile):
    path = root_path +  "/{pair}/{condition}/{node}".format(pair=pair,condition=condition,node=node)
    if sp == "dum":
        sister_sp = "ten"
    elif sp == "ten" or sp == "tent":
        sister_sp = "dum"
    elif sp == "sar":
        sister_sp = "bic"
    elif sp == "bic":
        sister_sp = "sar"
    gwf.target_from_template(
    name = "scaled_{quantile}_timeline_{node}".format(quantile=quantile,node=node),
    template = node_timeline_scaled(path,log_path,script_path,sp,sister_sp,node,quantile)
    )

def cal_scaled_terminal(gwf,root_path,log_path,script_path,sp,node,pair,condition,quantile):
    path = root_path +  "/{pair}/{condition}/{node}".format(pair=pair,condition=condition,node=node)
    if sp == "dum":
        sister_sp = "ten"
    elif sp == "ten" or sp == "tent":
        sister_sp = "dum"
    elif sp == "sar":
        sister_sp = "bic"
    elif sp == "bic":
        sister_sp = "sar"
    gwf.target_from_template(
    name = "scaled_{quantile}_terminal_{node}".format(quantile=quantile,node=node),
    template = node_terminal_scaled(path,log_path,script_path,sp,sister_sp,node,quantile)
    )

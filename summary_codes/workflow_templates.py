from gwf import *

def retrieve_independent_nodes(path,log_path,sp,node):
    inputs = [log_path+"/{sp}_{node}_blast_root.DONE".format(sp=sp,node=node)]
    outputs = [log_path + "/summary/{sp}_{node}_sum_prepare.DONE".format(sp=sp,node=node)]
    options = {
               'cores': 1,
               'memory': '1g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    mkdir -p {path}
    cd {path}
    conda activate biopython
    cp {log_path}/../../phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.fa {node}_aln.fa
    cp {log_path}/../../phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.dist {node}_aln.dist
    cp {log_path}/../../phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_rooted.newick {node}_rooted.newick
    cp {log_path}/../../phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_root.tsv {node}_root.tsv 
    echo done > {log}
    """.format(path=path,log_path=log_path,sp=sp,node=node,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def node_timeline_scaled(path,log_path,script_path,sp,sister_sp,node,quantile):
    inputs = [log_path + "/summary/{sp}_{node}_sum_prepare.DONE".format(sp=sp,node=node)]
    outputs = [log_path + "/summary/{sp}_{node}_scaled_{quantile}_timeline.DONE".format(sp=sp,node=node,quantile=quantile)] 
    options = {
               'cores': 1,
               'memory': '16g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    echo "getting scaled timeline"
    mkdir -p {path}
    cd {path}
    conda activate biopython
    Rscript {script_path}/rooted_node_vector_uniformed.R {node}_rooted.newick {node} {sp} {sister_sp} {quantile}
    echo done > {log}
    """.format(path=path,script_path = script_path, log = outputs[0],node=node,sp=sp,sister_sp=sister_sp,quantile = quantile) 
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def node_terminal_scaled(path,log_path,script_path,sp,sister_sp,node,quantile):
    inputs = [log_path + "/summary/{sp}_{node}_sum_prepare.DONE".format(sp=sp,node=node)]
    outputs = [log_path + "/summary/{sp}_{node}_scaled_{quantile}_terminal.DONE".format(sp=sp,node=node,quantile=quantile)] 
    options = {
               'cores': 1,
               'memory': '16g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    echo "getting scaled timeline"
    mkdir -p {path}
    cd {path}
    conda activate biopython
    Rscript {script_path}/rooted_terminal_time_uniformed.R {node}_rooted.newick {node} {sp} {sister_sp} {quantile}
    echo done > {log}
    """.format(path=path,script_path = script_path, log = outputs[0],node=node,sp=sp,sister_sp=sister_sp,quantile = quantile) 
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

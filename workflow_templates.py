from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/logs/orf"
def blast_makedb(path,fasta,name):
    inputs = [fasta]
    outputs = [LOG_PATH+"/{name}_blastdb.DONE".format(name=name)]
    options = {
               'cores': 4,
               'memory': '16g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate blast
    echo jobinfo $SLURM_JOBID
    echo "blast sv fasta files"
    mkdir -p {path}
    cd {path}
    ln -s {fasta} {name}.fa
    makeblastdb -in {name}.fa -dbtype nucl
    echo done > {log}
    """.format(path=path,fasta=fasta,name = name,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def blast_search(path,db,seed,outname):
    inputs = [LOG_PATH+"/{db}_blastdb.DONE".format(db=db),
            seed]
    outputs = [LOG_PATH+"/{outname}_blast.DONE".format(outname = outname)]
    options = {
               'cores':16,
               'memory': '16g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate blast
    echo jobinfo $SLURM_JOBID
    echo "blast sv fasta files"
    mkdir -p {path}
    cd {path}
    blastn -db /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/{db}.fa -num_threads 16 -query {seed} -out {outname}_blast.txt -outfmt "7 qseqid qstart qend sstart send pident length sseqid"
    echo done > {log}
    """.format(path=path,db=db,seed=seed,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def blast_search_node(path,db,seed,node,outname):
    inputs = [LOG_PATH+"/{db}_blastdb.DONE".format(db=db),
             LOG_PATH+"/build_seed_{node}.DONE".format(node=node)
             ]
    outputs = [LOG_PATH+"/{outname}_blast.DONE".format(outname = outname)]
    options = {
               'cores':16,
               'memory': '32g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate blast
    echo jobinfo $SLURM_JOBID
    echo "blast sv fasta files"
    mkdir -p {path}
    cd {path}
    blastn -db /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/blast/{db}.fa -num_threads 16 -query {seed} -out {outname}_blast.txt -outfmt "7 qseqid qstart qend sstart send pident length sseqid"
    echo done > {log}
    """.format(path=path,db=db,seed=seed,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def graph2uniq_node_seq(path,rds,social,subsocial):
    inputs = [rds]
    outputs = [LOG_PATH+"/{social}_{subsocial}_uniq_nodes_summary.DONE".format(social=social,subsocial=subsocial)]
    options = {
               'cores': 1,
               'memory': '4g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate Rseq
    echo jobinfo $SLURM_JOBID
    echo "blast sv fasta files"
    mkdir -p {path}
    cd {path}
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/graph2uniq_node_seq_id.R {rds} {social} {subsocial} 
    echo done > {log}
    """.format(path=path,rds=rds,social=social,subsocial=subsocial,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def make_node_txt(path,social,subsocial,seq_summary):
    inputs = [LOG_PATH+"/{social}_{subsocial}_uniq_nodes_summary.DONE".format(social=social,subsocial=subsocial)]
    outputs = [LOG_PATH+"/{social}_{subsocial}_per_node_summary.DONE".format(social=social,subsocial=subsocial)]
    options = {
               'cores': 1,
               'memory': '1g',
               'walltime':"1:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate Rseq
    echo jobinfo $SLURM_JOBID
    echo "blast sv fasta files"
    mkdir -p {path}
    cd {path}
    bash /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/process_seq_summary.sh {seq_summary}
    echo done > {log}
    """.format(path=path,seq_summary=seq_summary,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def unique_node_seq(path,social,subsocial):
    inputs = [LOG_PATH+"/{social}_{subsocial}_uniq_nodes_summary.DONE".format(social=social,subsocial=subsocial)]
    outputs = [LOG_PATH+"/{social}_{subsocial}_uniq_nodes_2fa.DONE".format(social=social,subsocial = subsocial)]
    options = {
               'cores': 1,
               'memory': '4g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    echo "create fasta files"
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/retrieve_rename_fasta.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/{social}_{subsocial}/{social}_{subsocial}_uniq_node_seq.tsv /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/pair_syri/{social}_{subsocial}/wgs/{social}_{subsocial}_wgs_comp2.fa /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/{social}_{subsocial}/{social}_{subsocial}_uniq_node_seq.fa
    echo done > {log}
    """.format(path=path,social=social,subsocial=subsocial,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def build_node_seed(path,pair,sp,node):
    inputs = ["/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/{pair}/{sp}/{node}_seq.txt".format(pair=pair,sp=sp,node=node)]
    outputs = [LOG_PATH+"/build_seed_{node}.DONE".format(node=node)]
    options = {
               'cores': 1,
               'memory': '16g',
               'walltime':"8:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    echo "create fasta files"
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/retrieve_fasta.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/{pair}/{pair}_uniq_node_seq.fa /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_0_seqs/{pair}/{sp}/{node}_seq.txt {node}_seq.fa
    echo done > {log}
    """.format(path=path,pair=pair,sp=sp,node=node,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def blast_compose(path,sp,node):
    inputs = [LOG_PATH + "/{sp}_{node}_blast.DONE".format(sp=sp,node=node)]
    outputs = [LOG_PATH + "/{sp}_{node}_blast_composed.DONE".format(sp=sp,node=node)]
    options = {
               'cores': 1,
               'memory': '64g',
               'walltime':"4:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo jobinfo $SLURM_JOBID
    echo "post processing blast hits"
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/node_blast_filter.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast.txt /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dup.txt
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/node_blast_filter_merge.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dup.txt /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass.txt 
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/blasthit2fa.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass.txt /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter3_dnm/data/genome/Ste_ncbi_chromosome.fa /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass.fa
    echo done > {log}
    """.format(path=path,sp=sp,node=node,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def align_blast(path,sp,node):
    inputs = [LOG_PATH + "/{sp}_{node}_blast_composed.DONE".format(sp=sp,node=node)]
    outputs = [LOG_PATH + "/{sp}_{node}_blast_aligned.DONE".format(sp=sp,node=node)]
    options = {
               'cores': 18,
               'memory': '200g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    echo "align blast hits"
    mkdir -p {path}
    cd {path}
    conda activate biopython
    echo "remove duplicated sequences"
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/fa_dedup.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass.fa /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup.fa 
    conda activate mafft
    echo "start alignment"
    mafft --thread 18 --threadtb 5 --threadit 0 --reorder --adjustdirection --auto /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup.fa > /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln.fa
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/aln_gap_filter.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln.fa /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.fa
    echo done > {log}
    """.format(path=path,sp=sp,node=node,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def build_aln_distance(path,sp,node):
    inputs = [LOG_PATH + "/{sp}_{node}_blast_aligned.DONE".format(sp=sp,node=node)]
    outputs = [LOG_PATH + "/{sp}_{node}_blast_dist.DONE".format(sp=sp,node=node)]
    options = {
               'cores': 1,
               'memory': '48g',
               'walltime':"48:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    echo "align blast hits to pairwise distance"
    mkdir -p {path}
    cd {path}
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/fa2dist_refine.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.fa 5 0 /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.dist
    echo done > {log}
    """.format(path=path,sp=sp,node=node,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def root_node(path,sp,node):
    inputs = [LOG_PATH + "/{sp}_{node}_blast_dist.DONE".format(sp=sp,node=node)]
    outputs = [LOG_PATH + "/{sp}_{node}_blast_root.DONE".format(sp=sp,node=node)]
    options = {
               'cores': 1,
               'memory': '32g',
               'walltime':"8:00:00",
               'account':"spider2"
               }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    echo "setting root information of the node"
    mkdir -p {path}
    cd {path}
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/root_set.py /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_1_aln/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.fa {sp} {sp}_{node}
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/dist2newick.R /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_blast_pass_dedup_aln_filtered.dist /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_root.tsv 
    echo done > {log}
    """.format(path=path,sp=sp,node=node,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def node_timeline(path,sp,sister_sp,node):
    inputs = [LOG_PATH + "/{sp}_{node}_blast_root.DONE".format(sp=sp,node=node)]
    outputs = [LOG_PATH + "/{sp}_{node}_blast_timeline.DONE".format(sp=sp,node=node)]
    options = {
               'cores': 1,
               'memory': '16g',
               'walltime':"2:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    echo "setting root information of the node"
    mkdir -p {path}
    cd {path}
    conda activate biopython
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/rooted_node_vector.R /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/step/phylo/step_2_phylo/{sp}/{node}/{sp}_{node}_rooted.newick {node} {sp} {sister_sp} 
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/chapter2_te/script/orf_finder/timeline_viz.R {node} {sp} {sister_sp}
    echo done > {log}
    """.format(path=path,sp=sp,node=node,log=outputs[0],sister_sp=sister_sp)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def raxml_tree(path,aln,outname):
    inputs = [aln]
    outputs = [LOG_PATH+"/raxml_{outname}.DONE".format(outname=outname)]
    options = {
               'cores': 16,
               'memory': '32g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate raxml
    echo jobinfo $SLURM_JOBID
    echo "raxml build ML tree with boostrapping"
    mkdir -p {path}
    cd {path}
    # Basic RAxML command for ML tree construction with bootstrap
    raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -p 12345 -x 12345 -T 16 -N 100 -s {aln} -n {outname} --print-identical-sequences
    # raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 12345 -T 16 -f b -t RAxML_bestTree.{outname} -z RAxML_bootstrap.{outname} -n {outname}_final --print-identical-sequences
    echo done > {log}
    """.format(path=path,aln = aln, outname = outname,log =outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

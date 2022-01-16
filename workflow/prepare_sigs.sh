snakemake -s prepare_sigs.smk --cluster "sbatch -A {cluster.account} -t {cluster.time} -J {cluster.jname} -c {cluster.ntasks-per-node} -p {cluster.partition} -N {cluster.nodes} --mem={cluster.mem}" --cluster-config cluster_config.yml --jobs 15
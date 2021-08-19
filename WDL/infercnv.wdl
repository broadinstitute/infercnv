version 1.0

workflow infercnv {
    input {
        File raw_counts_matrix # the matrix of genes (rows) vs. cells (columns)
        File gene_order_file # data file containing the positions of each gene along each chromosome in the genome
        File annotations_file # a description of the cells, indicating the cell type classifications.
        String additional_args = ""
        Int cpu = 1
        String memory = "12G"
        String docker = "trinityctat/infercnv:1.7.1-S"
        Int preemptible = 2
        Int extra_disk_space = 10
    }

    call run_infercnv {
        input:
            raw_counts_matrix=raw_counts_matrix,
            gene_order_file=gene_order_file,
            annotations_file=annotations_file,
            additional_args=additional_args,
            cpu=cpu,
            memory=memory,
            extra_disk_space=extra_disk_space,
            docker=docker,
            preemptible=preemptible
    }

    output {
        Array[File] infercnv_outputs = run_infercnv.infercnv_outputs
    }
}

task run_infercnv {
    input {
        File raw_counts_matrix
        File gene_order_file
        File annotations_file
        String memory
        Int cpu
        String docker
        Int preemptible
        String additional_args
        Int extra_disk_space
    }

    command <<<
        set -e

        mkdir infercnv

        inferCNV.R \
        --raw_counts_matrix ~{raw_counts_matrix} \
        --annotations_file ~{annotations_file} \
        --gene_order_file ~{gene_order_file} \
        --num_threads ~{cpu} \
        --out_dir infercnv \
        ~{additional_args} \
    >>>

    output {
        Array[File] infercnv_outputs = glob("infercnv/*")
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(raw_counts_matrix, "GB")*2 + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}


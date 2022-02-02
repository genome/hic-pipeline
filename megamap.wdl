version 1.0

import "./hic.wdl"

workflow megamap {
    meta {
        version: "1.10.0"
        caper_docker: "encodedcc/hic-pipeline:1.10.0"
        caper_singularity: "docker://encodedcc/hic-pipeline:1.10.0"
    }

    input {
        Array[File] bams
        Array[File] hic_files
        File? restriction_sites
        File chrom_sizes
        String assembly_name = "undefined"
        Int quality = 30
        Boolean intact = true
        File? phased_vcf
        Int? create_hic_num_cpus
        String delta_docker = "encodedcc/hic-pipeline:1.10.0_delta"
        String hiccups_docker = "encodedcc/hic-pipeline:1.10.0_hiccups"
    }

    String delta_models_path = if intact then "ultimate-models" else "beta-models"
    Array[Int] delta_resolutions = if intact then [5000, 2000, 1000] else [5000, 10000]
    Array[Int] create_hic_in_situ_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100]
    Array[Int] create_hic_intact_resolutions = [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000, 2000, 1000, 500, 200, 100, 50, 20, 10]
    Array[Int] create_hic_resolutions = if intact then create_hic_intact_resolutions else create_hic_in_situ_resolutions

    call hic.normalize_assembly_name as normalize_assembly_name { input:
        assembly_name = assembly_name
    }

    call hic.merge as merge { input:
        bams = bams,
    }

    call hic.bam_to_pre as bam_to_pre { input:
        bam = merge.bam,
        quality = quality,
    }

    call hic.create_accessibility_track as accessibility { input:
        pre = bam_to_pre.pre,
        chrom_sizes = chrom_sizes,
    }

    call merge_stats_from_hic_files { input:
        hic_files = hic_files,
    }

    if (normalize_assembly_name.assembly_is_supported) {
        call hic.create_hic as create_hic { input:
            pre = bam_to_pre.pre,
            pre_index = bam_to_pre.index,
            restriction_sites = restriction_sites,
            quality = quality,
            stats = merge_stats_from_hic_files.merged_stats,
            stats_hists = merge_stats_from_hic_files.merged_stats_hists,
            resolutions = create_hic_resolutions,
            assembly_name = normalize_assembly_name.normalized_assembly_name,
            num_cpus = create_hic_num_cpus,
        }
    }

    if (!normalize_assembly_name.assembly_is_supported) {
        call hic.create_hic as create_hic_with_chrom_sizes { input:
            pre = bam_to_pre.pre,
            pre_index = bam_to_pre.index,
            restriction_sites = restriction_sites,
            quality = quality,
            stats = merge_stats_from_hic_files.merged_stats,
            stats_hists = merge_stats_from_hic_files.merged_stats_hists,
            resolutions = create_hic_resolutions,
            assembly_name = assembly_name,
            num_cpus = create_hic_num_cpus,
            chrsz =  chrom_sizes,
        }
    }

    File unnormalized_hic_file = select_first([
        if (defined(create_hic.output_hic))
        then create_hic.output_hic
        else create_hic_with_chrom_sizes.output_hic
    ])

    call hic.add_norm as add_norm { input:
        hic = unnormalized_hic_file,
        quality = quality,
    }

    call hic.arrowhead as arrowhead { input:
        hic_file = add_norm.output_hic,
        quality = quality,
    }

    if (!intact) {
        call hic.hiccups { input:
            hic_file = add_norm.output_hic,
            quality = quality,
            docker = hiccups_docker,
        }
    }

    if (intact) {
        call hic.hiccups_2 { input:
            hic = add_norm.output_hic,
            quality = quality,
            docker = hiccups_docker,
        }

        call hic.localizer as localizer_intact { input:
            hic = add_norm.output_hic,
            loops = hiccups_2.merged_loops,
            quality = quality,
        }
    }

    call hic.create_eigenvector as create_eigenvector { input:
        hic_file = add_norm.output_hic,
        chrom_sizes = chrom_sizes,
        output_filename_suffix = "_" + quality,
    }

    call hic.create_eigenvector as create_eigenvector_10kb { input:
        hic_file = add_norm.output_hic,
        chrom_sizes = chrom_sizes,
        resolution = 10000,
        output_filename_suffix = "_" + quality,
    }

    call hic.delta as delta { input:
        hic = add_norm.output_hic,
        docker = delta_docker,
        resolutions = delta_resolutions,
        models_path = delta_models_path,
    }

    call hic.localizer as localizer_delta { input:
        hic = add_norm.output_hic,
        loops = delta.loops,
    }

    call hic.slice as slice_25kb { input:
        hic_file = add_norm.output_hic,
        resolution = 25000,
    }

    call hic.slice as slice_50kb { input:
        hic_file = add_norm.output_hic,
        resolution = 50000,
    }

    call hic.slice as slice_100kb { input:
        hic_file = add_norm.output_hic,
        resolution = 100000,
    }
}


task merge_stats_from_hic_files {
    input {
        Array[File] hic_files
        Int quality = 30
    }

    command <<<
        set -euo pipefail
        java \
            -jar \
            /opt/merge-stats.jar \
            ~{"inter_" + quality} \
            ~{sep=" " hic_files}
    >>>

    output {
        File merged_stats = "inter_~{quality}.txt"
        File merged_stats_hists = "inter_~{quality}_hists.m"
    }
}

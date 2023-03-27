from pathlib import Path

configfile: 'config.yml'
bam_folder = Path(config["BAMS_DIR"])

rule all:
    input:
        sample_manifesto = config["BAM_SAMPLES_JSON"]


def aggregate_downloaded_bams(wcards):
    # Halt workflow until BAMs are downloaded
    checkpoints.download_samples.get(**wcards)
    bam_files = [f for f in bam_folder.glob("*.ptrim.bam")]
    return bam_files


rule make_samples_list:
    input: aggregate_downloaded_bams
    output: config["BAM_SAMPLES_JSON"]
    conda: "envs/nywws_env.yml"
    shell:
        "python scripts/make_bam_json_samples.py {input} ; "


checkpoint download_samples:
    output: touch(bam_folder/".bams-downloaded")
    params:
        out = bam_folder
    shell:
        "rclone copy --progress gcs:su_nywws_test_bucket/nymc {params.out} --include \"/202*.ptrim.bam\" ; "
        "rclone copy --progress gcs:su_nywws_test_bucket/rochester/nywws {params.out} --include \"/202*.ptrim.bam\" ; "
        "rclone copy --progress gcs:su_nywws_test_bucket/suny_upstate {params.out} --include \"/202*.ptrim.bam\" ; "
        "rclone copy --progress gcs:su_nywws_test_bucket/wadsworth {params.out} --include \"/202*.ptrim.bam\" ; "

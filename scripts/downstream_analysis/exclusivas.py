#!/usr/bin/env python3

import os
import subprocess
import yaml
import shutil
from datetime import datetime


def run_cmd(cmd):
    print(f"[CMD] {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def load_samples(sample_list_path):
    samples = []
    with open(sample_list_path) as f:
        for line in f:
            if line.strip():
                name, path = line.strip().split()
                samples.append({"name": name, "vcf": path})
    return samples


def filter_vcf(sample, cfg, variant_type):
    outdir = cfg["paths"]["filtered_dir"]
    os.makedirs(outdir, exist_ok=True)

    min_dp = cfg["global_arguments"]["min_dp"]

    # tipo de variante
    if variant_type == "snps":
        type_flag = "-v snps"
    elif variant_type == "indels":
        type_flag = "-v indels"
    else:
        type_flag = ""

    output = f"{outdir}/{sample['name']}.{variant_type}.filtered.vcf.gz"
    sample["filtered"] = output

    # filtro PASS + tipo + DP mínimo
    cmd = (
        f"bcftools view -f PASS {type_flag} "
        f"-i 'FORMAT/DP>={min_dp}' {sample['vcf']} "
        f"| bcftools sort -Oz -o {output}"
    )
    run_cmd(cmd)

    run_cmd(f"bcftools index {output}")

    return sample


def compute_exclusive(samples, cfg, variant_type):
    outdir = cfg["paths"]["exclusive_dir"]
    os.makedirs(outdir, exist_ok=True)

    for s in samples:
        others = [x for x in samples if x["name"] != s["name"]]
        others_inputs = " ".join([o["filtered"] for o in others])

        tmp_dir = f"{outdir}/{s['name']}_{variant_type}_tmp"
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

        cmd = f"bcftools isec -C {s['filtered']} {others_inputs} -p {tmp_dir}"
        run_cmd(cmd)

        exclusive_vcf = f"{outdir}/{s['name']}.{variant_type}.exclusive.vcf.gz"
        s["exclusive"] = exclusive_vcf

        tmp_vcf = f"{tmp_dir}/0000.vcf"
        run_cmd(f"bcftools view {tmp_vcf} -Oz -o {exclusive_vcf}")
        run_cmd(f"bcftools index {exclusive_vcf}")

        shutil.rmtree(tmp_dir)

    return samples


def run_annovar(input_vcf, output_prefix, cfg):
    annovar_script = os.path.join(
        cfg["global_arguments"]["annovar_path"], "table_annovar.pl"
    )
    db = cfg["global_arguments"]["annovar_reference"]

    outdir = os.path.dirname(output_prefix)
    os.makedirs(outdir, exist_ok=True)

    cmd = (
        f"perl {annovar_script} {input_vcf} {db} "
        f"-buildver hg38 -out {output_prefix} "
        f"-remove -protocol refGene,avsnp150 "
        f"-operation g,f -nastring . -vcfinput"
    )
    run_cmd(cmd)


def main():
    print(f"Iniciando pipeline às {datetime.now().ctime()}\n")

    with open("config.yaml") as f:
        cfg = yaml.safe_load(f)

    variant_type = cfg["global_arguments"]["variant_type"]
    samples = load_samples(cfg["samples"]["sample_list"])

    print("=== FILTRANDO (PASS + tipo + DP mínimo) ===")
    samples = [filter_vcf(s, cfg, variant_type) for s in samples]

    print("\n=== EXCLUSIVAS ===")
    samples = compute_exclusive(samples, cfg, variant_type)

    print("\n=== ANNOTAÇÃO ANNOVAR (TOTAL + EXCLUSIVAS) ===")

    ann_total_dir = cfg["paths"]["annovar_total"]
    ann_excl_dir = cfg["paths"]["annovar_dir"]

    for s in samples:
        # Total PASS anotado
        run_annovar(
            s["filtered"],
            f"{ann_total_dir}/{s['name']}.{variant_type}.total",
            cfg,
        )

        # Exclusivas
        run_annovar(
            s["exclusive"],
            f"{ann_excl_dir}/{s['name']}.{variant_type}.exclusive",
            cfg,
        )

    print(f"\nPipeline finalizado às {datetime.now().ctime()}\n")


if __name__ == "__main__":
    main()


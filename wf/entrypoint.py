import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Annotated, List, Optional

import requests
from flytekit.core.annotation import FlyteAnnotation
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata

sys.stdout.reconfigure(line_buffering=True)


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    antibody: Optional[str] = None
    control: Optional[str] = None


class Reference_Type(Enum):
    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    run_name: Annotated[
        str,
        FlyteAnnotation(
            {
                "rules": [
                    {
                        "regex": r"^[a-zA-Z0-9_-]+$",
                        "message": "ID name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
                    }
                ],
            }
        ),
    ],
    pvc_name: str,
    input: List[SampleSheet],
    seq_center: Optional[str],
    read_length: Optional[int],
    genome_source: str,
    latch_genome: Reference_Type,
    outdir: LatchOutputDir,
    email: Optional[str],
    multiqc_title: Optional[str],
    genome: Optional[str],
    fasta: Optional[LatchFile],
    gtf: Optional[LatchFile],
    gff: Optional[LatchFile],
    bwa_index: Optional[str],
    bowtie2_index: Optional[str],
    chromap_index: Optional[str],
    star_index: Optional[str],
    gene_bed: Optional[LatchFile],
    macs_gsize: Optional[float],
    blacklist: Optional[str],
    save_reference: bool,
    clip_r1: Optional[int],
    clip_r2: Optional[int],
    three_prime_clip_r1: Optional[int],
    three_prime_clip_r2: Optional[int],
    trim_nextseq: Optional[int],
    skip_trimming: bool,
    save_trimmed: bool,
    keep_dups: bool,
    keep_multi_map: bool,
    bwa_min_score: Optional[int],
    save_align_intermeds: bool,
    save_unaligned: bool,
    narrow_peak: bool,
    macs_fdr: Optional[float],
    macs_pvalue: Optional[float],
    save_macs_pileup: bool,
    skip_peak_qc: bool,
    skip_peak_annotation: bool,
    skip_consensus_peaks: bool,
    skip_fastqc: bool,
    skip_picard_metrics: bool,
    skip_preseq: bool,
    skip_plot_profile: bool,
    skip_plot_fingerprint: bool,
    skip_spp: bool,
    skip_deseq2_qc: bool,
    skip_igv: bool,
    skip_multiqc: bool,
    skip_qc: bool,
    fragment_size: int,
    aligner: Optional[str],
    broad_cutoff: Optional[float],
    min_reps_consensus: Optional[int],
    deseq2_vst: bool,
) -> None:
    shared_dir = Path("/nf-workdir")
    rename_current_execution(str(run_name))

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    profile_list = ["docker", "test"]

    if len(profile_list) == 0:
        profile_list.append("standard")

    profiles = ",".join(profile_list)

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        profiles,
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input),
        *get_flag("fragment_size", fragment_size),
        *get_flag("seq_center", seq_center),
        *get_flag("read_length", read_length),
        *get_flag("outdir", LatchOutputDir(f"{outdir.remote_path}/{run_name}")),
        *get_flag("email", email),
        *get_flag("multiqc_title", multiqc_title),
        *get_flag("genome", genome),
        *get_flag("fasta", fasta),
        *get_flag("gtf", gtf),
        *get_flag("gff", gff),
        *get_flag("bwa_index", bwa_index),
        *get_flag("bowtie2_index", bowtie2_index),
        *get_flag("chromap_index", chromap_index),
        *get_flag("star_index", star_index),
        *get_flag("gene_bed", gene_bed),
        *get_flag("macs_gsize", macs_gsize),
        *get_flag("blacklist", blacklist),
        *get_flag("save_reference", save_reference),
        *get_flag("clip_r1", clip_r1),
        *get_flag("clip_r2", clip_r2),
        *get_flag("three_prime_clip_r1", three_prime_clip_r1),
        *get_flag("three_prime_clip_r2", three_prime_clip_r2),
        *get_flag("trim_nextseq", trim_nextseq),
        *get_flag("skip_trimming", skip_trimming),
        *get_flag("save_trimmed", save_trimmed),
        *get_flag("aligner", aligner),
        *get_flag("keep_dups", keep_dups),
        *get_flag("keep_multi_map", keep_multi_map),
        *get_flag("bwa_min_score", bwa_min_score),
        *get_flag("save_align_intermeds", save_align_intermeds),
        *get_flag("save_unaligned", save_unaligned),
        *get_flag("narrow_peak", narrow_peak),
        *get_flag("broad_cutoff", broad_cutoff),
        *get_flag("macs_fdr", macs_fdr),
        *get_flag("macs_pvalue", macs_pvalue),
        *get_flag("min_reps_consensus", min_reps_consensus),
        *get_flag("save_macs_pileup", save_macs_pileup),
        *get_flag("skip_peak_qc", skip_peak_qc),
        *get_flag("skip_peak_annotation", skip_peak_annotation),
        *get_flag("skip_consensus_peaks", skip_consensus_peaks),
        *get_flag("skip_fastqc", skip_fastqc),
        *get_flag("skip_picard_metrics", skip_picard_metrics),
        *get_flag("skip_preseq", skip_preseq),
        *get_flag("deseq2_vst", deseq2_vst),
        *get_flag("skip_plot_profile", skip_plot_profile),
        *get_flag("skip_plot_fingerprint", skip_plot_fingerprint),
        *get_flag("skip_spp", skip_spp),
        *get_flag("skip_deseq2_qc", skip_deseq2_qc),
        *get_flag("skip_igv", skip_igv),
        *get_flag("skip_multiqc", skip_multiqc),
        *get_flag("skip_qc", skip_qc),
    ]

    # Add genome-specific parameters if using Latch genome source
    if genome_source == "latch_genome_source":
        cmd += [
            "--fasta",
            f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.fna",
            "--gtf",
            f"s3://latch-public/nf-core/rnaseq/{latch_genome.name}/{latch_genome.name}.genomic.gtf",
        ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(
                    urljoins(
                        "latch:///your_log_dir/nf_nf_core_chipseq", name, "nextflow.log"
                    )
                )
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_chipseq(
    run_name: Annotated[
        str,
        FlyteAnnotation(
            {
                "rules": [
                    {
                        "regex": r"^[a-zA-Z0-9_-]+$",
                        "message": "ID name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
                    }
                ],
            }
        ),
    ],
    input: List[SampleSheet],
    seq_center: Optional[str],
    read_length: Optional[int],
    outdir: LatchOutputDir,
    genome_source: str,
    latch_genome: Reference_Type,
    email: Optional[str],
    multiqc_title: Optional[str],
    genome: Optional[str],
    fasta: Optional[LatchFile],
    gtf: Optional[LatchFile],
    gff: Optional[LatchFile],
    bwa_index: Optional[str],
    bowtie2_index: Optional[str],
    chromap_index: Optional[str],
    star_index: Optional[str],
    gene_bed: Optional[LatchFile],
    macs_gsize: Optional[float],
    blacklist: Optional[str],
    save_reference: bool,
    clip_r1: Optional[int],
    clip_r2: Optional[int],
    three_prime_clip_r1: Optional[int],
    three_prime_clip_r2: Optional[int],
    trim_nextseq: Optional[int],
    skip_trimming: bool,
    save_trimmed: bool,
    keep_dups: bool,
    keep_multi_map: bool,
    bwa_min_score: Optional[int],
    save_align_intermeds: bool,
    save_unaligned: bool,
    narrow_peak: bool,
    macs_fdr: Optional[float],
    macs_pvalue: Optional[float],
    save_macs_pileup: bool,
    skip_peak_qc: bool,
    skip_peak_annotation: bool,
    skip_consensus_peaks: bool,
    skip_fastqc: bool,
    skip_picard_metrics: bool,
    skip_preseq: bool,
    skip_plot_profile: bool,
    skip_plot_fingerprint: bool,
    skip_spp: bool,
    skip_deseq2_qc: bool,
    skip_igv: bool,
    skip_multiqc: bool,
    skip_qc: bool,
    fragment_size: int = 200,
    aligner: Optional[str] = "bwa",
    broad_cutoff: Optional[float] = 0.1,
    min_reps_consensus: Optional[int] = 1,
    deseq2_vst: bool = True,
) -> None:
    """
    nf-core/chipseq

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(
        run_name=run_name,
        pvc_name=pvc_name,
        input=input,
        fragment_size=fragment_size,
        seq_center=seq_center,
        read_length=read_length,
        outdir=outdir,
        genome_source=genome_source,
        latch_genome=latch_genome,
        email=email,
        multiqc_title=multiqc_title,
        genome=genome,
        fasta=fasta,
        gtf=gtf,
        gff=gff,
        bwa_index=bwa_index,
        bowtie2_index=bowtie2_index,
        chromap_index=chromap_index,
        star_index=star_index,
        gene_bed=gene_bed,
        macs_gsize=macs_gsize,
        blacklist=blacklist,
        save_reference=save_reference,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        trim_nextseq=trim_nextseq,
        skip_trimming=skip_trimming,
        save_trimmed=save_trimmed,
        aligner=aligner,
        keep_dups=keep_dups,
        keep_multi_map=keep_multi_map,
        bwa_min_score=bwa_min_score,
        save_align_intermeds=save_align_intermeds,
        save_unaligned=save_unaligned,
        narrow_peak=narrow_peak,
        broad_cutoff=broad_cutoff,
        macs_fdr=macs_fdr,
        macs_pvalue=macs_pvalue,
        min_reps_consensus=min_reps_consensus,
        save_macs_pileup=save_macs_pileup,
        skip_peak_qc=skip_peak_qc,
        skip_peak_annotation=skip_peak_annotation,
        skip_consensus_peaks=skip_consensus_peaks,
        skip_fastqc=skip_fastqc,
        skip_picard_metrics=skip_picard_metrics,
        skip_preseq=skip_preseq,
        deseq2_vst=deseq2_vst,
        skip_plot_profile=skip_plot_profile,
        skip_plot_fingerprint=skip_plot_fingerprint,
        skip_spp=skip_spp,
        skip_deseq2_qc=skip_deseq2_qc,
        skip_igv=skip_igv,
        skip_multiqc=skip_multiqc,
        skip_qc=skip_qc,
    )

from pathlib import Path
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.services.register.utils import import_module_by_path

from wf.entrypoint import Reference_Type, SampleSheet, initialize, nextflow_runtime

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)


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
    nfcore/chipseq is a bioinformatics analysis pipeline used for Chromatin ImmunopreciPitation sequencing (ChIP-seq) data.


    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">

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


LaunchPlan(
    nf_nf_core_chipseq,
    "test_data",
    {
        "run_name": "Test-1",
    },
)

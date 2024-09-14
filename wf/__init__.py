from pathlib import Path
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.services.register.utils import import_module_by_path

from wf.entrypoint import (
    Aligner,
    Reference_Type,
    SampleSheet,
    initialize,
    nextflow_runtime,
)

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_chipseq(
    run_name: str,
    input: List[SampleSheet],
    seq_center: Optional[str],
    read_length: Optional[int],
    genome_source: str,
    latch_genome: Reference_Type,
    email: Optional[str],
    multiqc_title: Optional[str],
    genome: Optional[str],
    fasta: Optional[LatchFile],
    gtf: Optional[LatchFile],
    gff: Optional[LatchFile],
    bwa_index: Optional[LatchDir],
    bowtie2_index: Optional[LatchDir],
    chromap_index: Optional[LatchFile],
    star_index: Optional[LatchDir],
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
    aligner: Aligner = Aligner.bwa,
    broad_cutoff: Optional[float] = 0.1,
    min_reps_consensus: Optional[int] = 1,
    deseq2_vst: bool = True,
    outdir: LatchOutputDir = LatchOutputDir("latch:///Chipseq"),
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

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.3240506-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.3240506)

    **nfcore/chipseq** is a bioinformatics analysis pipeline used for Chromatin ImmunopreciPitation sequencing (ChIP-seq) data.

    This workflow is hosted on Latch Workflows, using a native Nextflow integration, with a graphical interface for accessible analysis by scientists. There is also an integration with Latch Registry so that batched workflows can be launched from “graphical sample sheets” or tables associating raw sequencing files with metadata.

    ## Documentation

    The nf-core/chipseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/chipseq/usage) and [output](https://nf-co.re/chipseq/output).

    The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

    ## Online videos

    A short talk about the history, current status and functionality on offer in this pipeline was given by Jose Espinosa-Carrasco ([@joseespinosa](https://github.com/joseespinosa)) on [26th July 2022](https://nf-co.re/events/2022/bytesize-chipseq) as part of the nf-core/bytesize series.

    You can find numerous talks on the [nf-core events page](https://nf-co.re/events) from various topics including writing pipelines/modules in Nextflow DSL2, using nf-core tooling, running nf-core pipelines as well as more generic content like contributing to Github. Please check them out!

    ## Pipeline summary

    1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
    3. Choice of multiple aligners
       1.([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
       2.([`Chromap`](https://github.com/haowenz/chromap)). **For paired-end reads only working until mapping steps, see [here](https://github.com/nf-core/chipseq/issues/291)**
       3.([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
       4.([`STAR`](https://github.com/alexdobin/STAR))
    4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
       1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
       2. Filtering to remove:
          - reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
          - reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
          - reads that are not marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
          - reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
          - reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
          - reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
          - reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); _paired-end only_)
          - reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
          - reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
          - reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); _paired-end only_)
       3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
       4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
       5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
       6. Calculate genome-wide IP enrichment relative to control ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
       7. Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))
       8. Call broad/narrow peaks ([`MACS2`](https://github.com/macs3-project/MACS))
       9. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
       10. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
       11. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
       12. PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
    6. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
    7. Present QC for raw read, alignment, peak-calling and differential binding results ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))

    ## Credits

    These scripts were originally written by Chuan Wang ([@chuan-wang](https://github.com/chuan-wang)) and Phil Ewels ([@ewels](https://github.com/ewels)) for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden. The pipeline was re-implemented by Harshil Patel ([@drpatelh](https://github.com/drpatelh)) from [Seqera Labs, Spain](https://seqera.io/) and converted to Nextflow DSL2 by Jose Espinosa-Carrasco ([@JoseEspinosa](https://github.com/JoseEspinosa)) from [The Comparative Bioinformatics Group](https://www.crg.eu/en/cedric_notredame) at [The Centre for Genomic Regulation, Spain](https://www.crg.eu/).

    Many thanks to others who have helped out and contributed along the way too, including (but not limited to): [@apeltzer](https://github.com/apeltzer), [@bc2zb](https://github.com/bc2zb), [@crickbabs](https://github.com/crickbabs), [@drejom](https://github.com/drejom), [@houghtos](https://github.com/houghtos), [@KevinMenden](https://github.com/KevinMenden), [@mashehu](https://github.com/mashehu), [@pditommaso](https://github.com/pditommaso), [@Rotholandus](https://github.com/Rotholandus), [@sofiahaglund](https://github.com/sofiahaglund), [@tiagochst](https://github.com/tiagochst) and [@winni2k](https://github.com/winni2k).

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#chipseq` channel](https://nfcore.slack.com/channels/chipseq) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    If you use nf-core/chipseq for your analysis, please cite it using the following doi: [10.5281/zenodo.3240506](https://doi.org/10.5281/zenodo.3240506)

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize(run_name=run_name)
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
    "Test Data",
    {
        "input": [
            SampleSheet(
                sample="SPT5_T0_REP1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822153_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822153_2.fastq.gz"
                ),
                antibody="SPT5",
                control="SPT5_INPUT_REP1",
            ),
            SampleSheet(
                sample="SPT5_T0_REP2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822154_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822154_2.fastq.gz"
                ),
                antibody="SPT5",
                control="SPT5_INPUT_REP2",
            ),
            SampleSheet(
                sample="SPT5_T15_REP1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822157_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822157_2.fastq.gz"
                ),
                antibody="SPT5",
                control="SPT5_INPUT_REP1",
            ),
            SampleSheet(
                sample="SPT5_T15_REP2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822158_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR1822158_2.fastq.gz"
                ),
                antibody="SPT5",
                control="SPT5_INPUT_REP2",
            ),
            SampleSheet(
                sample="SPT5_INPUT_REP1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR5204809_Spt5-ChIP_Input1_SacCer_ChIP-Seq_ss100k_R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR5204809_Spt5-ChIP_Input1_SacCer_ChIP-Seq_ss100k_R2.fastq.gz"
                ),
                antibody=None,
                control=None,
            ),
            SampleSheet(
                sample="SPT5_INPUT_REP2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR5204810_Spt5-ChIP_Input2_SacCer_ChIP-Seq_ss100k_R1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/chipseq/test_data/SRR5204810_Spt5-ChIP_Input2_SacCer_ChIP-Seq_ss100k_R2.fastq.gz"
                ),
                antibody=None,
                control=None,
            ),
        ],
        "run_name": "Test_Run",
        "genome_source": "custom",
        "read_length": 50,
        "skip_preseq": True,
        "fasta": LatchFile("s3://latch-public/nf-core/chipseq/test_data/genome.fa"),
        "gtf": LatchFile("s3://latch-public/nf-core/chipseq/test_data/genes.gtf"),
    },
)

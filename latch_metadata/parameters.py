from dataclasses import dataclass
from enum import Enum
from typing import Annotated, List, Optional

from flytekit.core.annotation import FlyteAnnotation
from latch.types.directory import LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchAuthor,
    NextflowMetadata,
    NextflowParameter,
    NextflowRuntimeResources,
    Params,
    Section,
    Spoiler,
    Text,
)


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    antibody: Optional[str]
    control: Optional[str]


class Reference_Type(Enum):
    homo_sapiens = "Homo sapiens (RefSeq GRCh38.p14)"
    mus_musculus = "Mus musculus (RefSeq GRCm39)"
    rattus_norvegicus = "Rattus norvegicus (RefSeq GRCr8)"
    # drosophila_melanogaster = "Drosophila melanogaster (RefSeq Release_6_plus_ISO1_MT)"
    # rhesus_macaque = "Macaca mulatta (RefSeq rheMac10/Mmul_10)"
    # saccharomyces_cerevisiae = "Saccharomyces cerevisiae (RefSeq R64)"


flow = [
    Section(
        "Input",
        Params(
            "input",
            "fragment_size",
        ),
    ),
    Section(
        "Reference Genome",
        Fork(
            "genome_source",
            "",
            latch_genome_source=ForkBranch(
                "Latch Verified Reference Genome",
                Params(
                    "latch_genome",
                ),
            ),
            custom=ForkBranch(
                "Custom Reference Genome",
                Params(
                    "fasta",
                    "gtf",
                    "gff",
                ),
                Spoiler(
                    "Advanced Options",
                    Params(
                        "bwa_index",
                        "bowtie2_index",
                        "chromap_index",
                        "star_index",
                        "gene_bed",
                        "macs_gsize",
                        "blacklist",
                        "save_reference",
                    ),
                ),
            ),
        ),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "Advanced Options",
        Section(
            "Adapter Trimming Options",
            Params(
                "clip_r1",
                "clip_r2",
                "three_prime_clip_r1",
                "three_prime_clip_r2",
                "trim_nextseq",
                "skip_trimming",
                "save_trimmed",
            ),
        ),
        Section(
            "Alignment Options",
            Params(
                "aligner",
                "keep_dups",
                "keep_multi_map",
                "bwa_min_score",
                "save_align_intermeds",
                "save_unaligned",
            ),
        ),
        Section(
            "Peak Calling Options",
            Params(
                "read_length",
                "narrow_peak",
                "broad_cutoff",
                "macs_fdr",
                "macs_pvalue",
                "min_reps_consensus",
                "save_macs_pileup",
                "skip_peak_qc",
                "skip_peak_annotation",
                "skip_consensus_peaks",
            ),
        ),
        Section(
            "Advanced customization Options",
            Params("seq_center", "email", "multiqc_title"),
        ),
        Section(
            "Skip processes",
            Params(
                "skip_fastqc",
                "skip_picard_metrics",
                "skip_preseq",
                "deseq2_vst",
                "skip_plot_profile",
                "skip_plot_fingerprint",
                "skip_spp",
                "skip_deseq2_qc",
                "skip_igv",
                "skip_multiqc",
                "skip_qc",
            ),
        ),
    ),
]


generated_parameters = {
    "run_name": NextflowParameter(
        type=str,
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
    ),
    "input": NextflowParameter(
        type=List[SampleSheet],
        display_name="Sample Sheet",
        samplesheet=True,
        samplesheet_type="csv",
        description="Samplesheet containing information about the samples in the experiment.",
    ),
    "fragment_size": NextflowParameter(
        type=int,
        display_name="Fragment Size",
        default=200,
        section_title=None,
        description="Estimated fragment size used to extend single-end reads.",
    ),
    "seq_center": NextflowParameter(
        type=Optional[str],
        display_name="Sequencing Center",
        default=None,
        section_title=None,
        description="Sequencing center information to be added to read group of BAM files.",
    ),
    "read_length": NextflowParameter(
        type=Optional[int],
        display_name="Read Length",
        default=None,
        section_title=None,
        description="Read length used to calculate MACS2 genome size for peak calling if `--macs_gsize` isn't provided.",
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        display_name="Output Directory",
        default=None,
        section_title=None,
        description="The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.",
    ),
    "genome_source": NextflowParameter(
        type=str,
        display_name="Reference Genome",
        description="Choose Reference Genome",
    ),
    "latch_genome": NextflowParameter(
        type=Reference_Type,
        display_name="Latch Verfied Reference Genome",
        description="Name of Latch Verfied Reference Genome.",
        default=Reference_Type.homo_sapiens,
    ),
    "email": NextflowParameter(
        type=Optional[str],
        display_name="Email Address",
        default=None,
        section_title=None,
        description="Email address for completion summary.",
    ),
    "multiqc_title": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Report Title",
        default=None,
        section_title=None,
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "genome": NextflowParameter(
        type=Optional[str],
        display_name="Reference Genome",
        default=None,
        description="Name of iGenomes reference.",
    ),
    "fasta": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Genome FASTA File",
        default=None,
        section_title=None,
        description="Path to FASTA genome file.",
    ),
    "gtf": NextflowParameter(
        type=Optional[LatchFile],
        display_name="GTF Annotation File",
        default=None,
        section_title=None,
        description="Path to GTF annotation file.",
    ),
    "gff": NextflowParameter(
        type=Optional[LatchFile],
        display_name="GFF3 Annotation File",
        default=None,
        section_title=None,
        description="Path to GFF3 annotation file.",
    ),
    "bwa_index": NextflowParameter(
        type=Optional[str],
        display_name="BWA Index",
        default=None,
        section_title=None,
        description="Path to directory or tar.gz archive for pre-built BWA index.",
    ),
    "bowtie2_index": NextflowParameter(
        type=Optional[str],
        display_name="Bowtie2 Index",
        default=None,
        section_title=None,
        description="Path to directory or tar.gz archive for pre-built Bowtie2 index.",
    ),
    "chromap_index": NextflowParameter(
        type=Optional[str],
        display_name="Chromap Index",
        default=None,
        section_title=None,
        description="Path to directory or tar.gz archive for pre-built Chromap index.",
    ),
    "star_index": NextflowParameter(
        type=Optional[str],
        display_name="STAR Index",
        default=None,
        section_title=None,
        description="Path to directory or tar.gz archive for pre-built STAR index.",
    ),
    "gene_bed": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Gene BED File",
        default=None,
        section_title=None,
        description="Path to BED file containing gene intervals. This will be created from the GTF file if not specified.",
    ),
    "macs_gsize": NextflowParameter(
        type=Optional[float],
        display_name="MACS2 Genome Size",
        default=None,
        section_title=None,
        description="Effective genome size parameter required by MACS2.",
    ),
    "blacklist": NextflowParameter(
        type=Optional[str],
        display_name="Blacklist Regions",
        default=None,
        section_title=None,
        description="Path to blacklist regions in BED format, used for filtering alignments.",
    ),
    "save_reference": NextflowParameter(
        type=bool,
        display_name="Save Reference",
        default=None,
        section_title=None,
        description="If generated by the pipeline save the BWA index in the results directory.",
    ),
    "clip_r1": NextflowParameter(
        type=Optional[int],
        display_name="Clip R1",
        default=None,
        description="Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).",
    ),
    "clip_r2": NextflowParameter(
        type=Optional[int],
        display_name="Clip R2",
        default=None,
        section_title=None,
        description="Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).",
    ),
    "three_prime_clip_r1": NextflowParameter(
        type=Optional[int],
        display_name="3' Clip R1",
        default=None,
        section_title=None,
        description="Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed.",
    ),
    "three_prime_clip_r2": NextflowParameter(
        type=Optional[int],
        display_name="3' Clip R2",
        default=None,
        section_title=None,
        description="Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed.",
    ),
    "trim_nextseq": NextflowParameter(
        type=Optional[int],
        display_name="Trim NextSeq",
        default=None,
        section_title=None,
        description="Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails.",
    ),
    "skip_trimming": NextflowParameter(
        type=bool,
        display_name="Skip Trimming",
        default=None,
        section_title=None,
        description="Skip the adapter trimming step.",
    ),
    "save_trimmed": NextflowParameter(
        type=bool,
        display_name="Save Trimmed Reads",
        default=None,
        section_title=None,
        description="Save the trimmed FastQ files in the results directory.",
    ),
    "aligner": NextflowParameter(
        type=Optional[str],
        display_name="Aligner",
        default="bwa",
        description="Specifies the alignment algorithm to use - available options are 'bwa', 'bowtie2' and 'star'.",
    ),
    "keep_dups": NextflowParameter(
        type=bool,
        display_name="Keep Duplicates",
        default=None,
        section_title=None,
        description="Duplicate reads are not filtered from alignments.",
    ),
    "keep_multi_map": NextflowParameter(
        type=bool,
        display_name="Keep Multi-mapped Reads",
        default=None,
        section_title=None,
        description="Reads mapping to multiple locations are not filtered from alignments.",
    ),
    "bwa_min_score": NextflowParameter(
        type=Optional[int],
        display_name="BWA Minimum Score",
        default=None,
        section_title=None,
        description="Don't output BWA MEM alignments with score lower than this parameter.",
    ),
    "save_align_intermeds": NextflowParameter(
        type=bool,
        display_name="Save Alignment Intermediates",
        default=None,
        section_title=None,
        description="Save the intermediate BAM files from the alignment step.",
    ),
    "save_unaligned": NextflowParameter(
        type=bool,
        display_name="Save Unaligned Reads",
        default=None,
        section_title=None,
        description="Where possible, save unaligned reads from either STAR, HISAT2 or Salmon to the results directory.",
    ),
    "narrow_peak": NextflowParameter(
        type=bool,
        display_name="Narrow Peak Mode",
        default=None,
        description="Run MACS2 in narrowPeak mode.",
    ),
    "broad_cutoff": NextflowParameter(
        type=Optional[float],
        display_name="Broad Cutoff",
        default=0.1,
        section_title=None,
        description="Specifies broad cutoff value for MACS2. Only used when --narrow_peak isnt specified.",
    ),
    "macs_fdr": NextflowParameter(
        type=Optional[float],
        display_name="MACS2 FDR",
        default=None,
        section_title=None,
        description="Minimum FDR (q-value) cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive.",
    ),
    "macs_pvalue": NextflowParameter(
        type=Optional[float],
        display_name="MACS2 p-value",
        default=None,
        section_title=None,
        description="p-value cutoff for peak detection, --macs_fdr and --macs_pvalue are mutually exclusive. If --macs_pvalue cutoff is set, q-value will not be calculated and reported as -1 in the final .xls file.",
    ),
    "min_reps_consensus": NextflowParameter(
        type=Optional[int],
        display_name="Minimum Replicates for Consensus",
        default=1,
        section_title=None,
        description="Number of biological replicates required from a given condition for a peak to contribute to a consensus peak.",
    ),
    "save_macs_pileup": NextflowParameter(
        type=bool,
        display_name="Save MACS2 Pileup",
        default=None,
        section_title=None,
        description="Instruct MACS2 to create bedGraph files normalised to signal per million reads.",
    ),
    "skip_peak_qc": NextflowParameter(
        type=bool,
        display_name="Skip Peak QC",
        default=None,
        section_title=None,
        description="Skip MACS2 peak QC plot generation.",
    ),
    "skip_peak_annotation": NextflowParameter(
        type=bool,
        display_name="Skip Peak Annotation",
        default=None,
        section_title=None,
        description="Skip annotation of MACS2 and consensus peaks with HOMER.",
    ),
    "skip_consensus_peaks": NextflowParameter(
        type=bool,
        display_name="Skip Consensus Peaks",
        default=None,
        section_title=None,
        description="Skip consensus peak generation, annotation and counting.",
    ),
    "skip_fastqc": NextflowParameter(
        type=bool,
        display_name="Skip FastQC",
        default=None,
        description="Skip FastQC.",
    ),
    "skip_picard_metrics": NextflowParameter(
        type=bool,
        display_name="Skip Picard Metrics",
        default=None,
        section_title=None,
        description="Skip Picard CollectMultipleMetrics.",
    ),
    "skip_preseq": NextflowParameter(
        type=bool,
        display_name="Skip Preseq",
        default=None,
        section_title=None,
        description="Skip Preseq.",
    ),
    "deseq2_vst": NextflowParameter(
        type=bool,
        display_name="DESeq2 VST",
        default=True,
        section_title=None,
        description="Use vst transformation instead of rlog with DESeq2.",
    ),
    "skip_plot_profile": NextflowParameter(
        type=bool,
        display_name="Skip Plot Profile",
        default=None,
        section_title=None,
        description="Skip deepTools plotProfile.",
    ),
    "skip_plot_fingerprint": NextflowParameter(
        type=bool,
        display_name="Skip Plot Fingerprint",
        default=None,
        section_title=None,
        description="Skip deepTools plotFingerprint.",
    ),
    "skip_spp": NextflowParameter(
        type=bool,
        display_name="Skip SPP",
        default=None,
        section_title=None,
        description="Skip Phantompeakqualtools.",
    ),
    "skip_deseq2_qc": NextflowParameter(
        type=bool,
        display_name="Skip DESeq2 QC",
        default=None,
        section_title=None,
        description="Skip DESeq2 PCA and heatmap plotting.",
    ),
    "skip_igv": NextflowParameter(
        type=bool,
        display_name="Skip IGV",
        default=None,
        section_title=None,
        description="Skip IGV.",
    ),
    "skip_multiqc": NextflowParameter(
        type=bool,
        display_name="Skip MultiQC",
        default=None,
        section_title=None,
        description="Skip MultiQC.",
    ),
    "skip_qc": NextflowParameter(
        type=bool,
        display_name="Skip QC",
        default=None,
        section_title=None,
        description="Skip all QC steps except for MultiQC.",
    ),
}

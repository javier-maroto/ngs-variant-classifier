# See https://varlociraptor.github.io/docs/calling/ for reference

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --id ID            Set the ID
  --reference FILE   Reference file (e.g., hg19.fa)
  --bam FILE         Input BAM file
  --vcf FILE         Input VCF file with candidates
  --n_cores N        Number of cores to use
  --max_depth N      Maximum depth value (if lower than DP, Varlociraptor downsamples before estimating)
  --use_mapq         Use MAPQ in Varlociraptor estimations
  --tumor_cell_quantity FLOAT      Tumor cell quantity in the sample (ranges from 0 to 1)
  --output_folder    Path where the outputs (data, logs and results) will be saved. Default: output
  --scenario         FFPE scenario (default: "ffpe"). "ffpe" to use contamination, "ffpe_no_cont" to not use contamination (--tumor_cell_quantity ignored in this case)
  --help             Show this help and exit
EOF
}

use_mapq=false
scenario="ffpe"
tumor_cell_quantity=1.0
output_folder="output"
n_cores=1
max_depth=20000

while [[ $# -gt 0 ]]; do
  case $1 in
    --id)   id="$2"; shift 2 ;;
    --reference)   reference="$2"; shift 2 ;;
    --bam)  bam="$2"; shift 2 ;;
    --vcf)  vcf="$2"; shift 2 ;;
    --n_cores)  n_cores=$2; shift 2 ;;
    --max_depth)  max_depth=$2; shift 2 ;;
    --use_mapq)   use_mapq=true; shift 1 ;;
    --tumor_cell_quantity) tumor_cell_quantity="$2"; shift 2 ;;
    --output_folder) output_folder="$2"; shift 2 ;;
    --scenario)
      scenario="$2"
      if [[ "$scenario" != "ffpe" && "$scenario" != "ffpe_no_cont" ]]; then
        echo "Error: --scenario must be either 'ffpe' or 'ffpe_no_cont'"
        usage
        exit 1
      fi
      shift 2
      ;;
    --help)       usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

# Make folders
mkdir -p ${output_folder}/${id}/data/chunks
mkdir -p ${output_folder}/${id}/data/observations
mkdir -p ${output_folder}/${id}/logs/preprocess
mkdir -p ${output_folder}/${id}/results

# Create scenario file
if [[ "$scenario" == "ffpe" ]]; then
  normal_score=$(awk -v v="$tumor_cell_quantity" 'BEGIN{print 1-v}')
  sed -r 's/^(\s*)(fraction\s*:\s*$)/\1fraction: '"$normal_score"'/' scenarios/scenario_ffpe.yaml > scenarios/scenario.yaml
else
  cp scenarios/scenario_ffpe_no_cont.yaml scenarios/scenario.yaml
  tumor_cell_quantity=1.0
fi

# Fix invalid records
bcftools norm -c s -m- -f $reference $vcf -Ov -o ${output_folder}/${id}/data/candidates_all.vcf > /dev/null

# Filter records
bcftools view -e 'FMT/AO==0' ${output_folder}/${id}/data/candidates_all.vcf -o ${output_folder}/${id}/data/candidates.vcf > /dev/null

# Convert to tsv
(
  echo -e "seqnames\tposition\tref\talt\tQUAL\tFILTER\tAF\tAO\tDP\tFAO\tFDP"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%AF\t%AO\t%DP\t%FAO\t%FDP]\n' \
    ${output_folder}/${id}/data/candidates.vcf
) > ${output_folder}/${id}/data/candidates.tsv

# Estimate alignment properties. 
# It assumes 
#   1. the set of samples has been sequenced and prepared in exactly the same way.
#   2. that MAPQ is as accurate as possible, and that too large indels are encoded as softclips (BWA mem style)
#   3. the data is not paired end
#     If not 3: consider manually providing --alignment-properties, e.g. computed with `samtools stats`.
time varlociraptor estimate alignment-properties $reference --bams $bam \
> ${output_folder}/${id}/data/alignment-properties.json 2>> ${output_folder}/${id}/logs/alignment-properties.log

# Build chunk names
chunk_names=()
for i in $(seq 1 $n_cores); do
    chunk_names+=("${output_folder}/${id}/data/chunks/chunk${i}.bcf")
done

# Split candidates. It assumes no structural variants (sorting maybe needed in that case)
rbt vcf-split ${output_folder}/${id}/data/candidates.vcf "${chunk_names[@]}"

# If Ion Torrent amplicons are unique and multi-mapping is rare, not using BAM in BWA style format might be fine
# It assumes 
#   1. the candidates are going to be the same across samples (by joint calling with a variant calling or merge candidates from different samples)
#   2. there are homopolymer sequencing errors such as Oxford Nanopore Technologies (ONT) reads
#     - If not 2: do not use argument --pairhmm-mode homopolymer
#   3. estimated AF is important, we want to minimize single read errors, and there are discrete allele frequencies of interest (e.g., 0.5, 1.0)
#     - If not 3: use argument --pairhmm-mode fast , for much faster computation (quadratic -> linear)
#   4. atomic candidate variants to deactivate realignment for SNVs and MNVs in varlociraptor
#     - If not 4: remove flag --atomic-candidate-variants if given candidate variants are merged into small haplotypes (given as MNVs or more complex replacements like provided e.g. by freebayes)
#   5. that it is processing amplicon data, where indels do impact the observed insert size
#     - If not 5: add flag --omit-insert-size
#   6. that MAPQ can be safely ignored
#     - If not 6: remove flag --omit-mapq-adjustment
if [[ "$use_mapq" == true ]]; then
  for f in $(seq 1 $n_cores); do
  (
  varlociraptor preprocess variants $reference --atomic-candidate-variants  \
  --pairhmm-mode homopolymer \
  --alignment-properties ${output_folder}/${id}/data/alignment-properties.json \
  --candidates ${output_folder}/${id}/data/chunks/chunk${f}.bcf \
  --bam $bam \
  --max-depth $max_depth \
  --output ${output_folder}/${id}/data/observations/chunk${f}.bcf 2>> ${output_folder}/${id}/logs/preprocess/chunk${f}.log
  bcftools index ${output_folder}/${id}/data/observations/chunk${f}.bcf
  ) &
  done
  wait
else
  for f in $(seq 1 $n_cores); do
  (
  varlociraptor preprocess variants $reference --atomic-candidate-variants --omit-mapq-adjustment \
  --pairhmm-mode homopolymer \
  --alignment-properties ${output_folder}/${id}/data/alignment-properties.json \
  --candidates ${output_folder}/${id}/data/chunks/chunk${f}.bcf \
  --bam $bam \
  --max-depth $max_depth \
  --output ${output_folder}/${id}/data/observations/chunk${f}.bcf 2>> ${output_folder}/${id}/logs/preprocess/chunk${f}.log
  bcftools index ${output_folder}/${id}/data/observations/chunk${f}.bcf
  ) &
  done
  wait
fi


# Merge results. Assumes there are no structural variants (option -a)
bcftools concat -a -Ob -o ${output_folder}/${id}/data/observations.bcf ${output_folder}/${id}/data/observations/chunk*.bcf
bcftools index ${output_folder}/${id}/data/observations.bcf


# For Ion Torrent amplicons, we use --omit-read-position-bias to avoid detecting a read position bias for those variants covered only by a single amplicon.
# It assumes:
#   1. that the variant caller found loci are correct (no de-novo sequencing, and MAPQ good enough)
#     - If not 1: remove flag --omit-alt-locus-bias
time varlociraptor call variants --omit-read-position-bias --omit-alt-locus-bias generic --scenario scenarios/scenario.yaml \
  --obs ffpetumor=${output_folder}/${id}/data/observations.bcf  > ${output_folder}/${id}/data/calls_filt.bcf 2>> ${output_folder}/${id}/logs/calling.log

bcftools view ${output_folder}/${id}/data/calls_filt.bcf -Ov -o ${output_folder}/${id}/data/calls_filt.vcf

(
  echo -e "seqnames\tposition\tref\talt\tQUAL\tFILTER\tAF\tDP\tSAOBS\tSROBS\tOOBS\tSB\tROB\tRPB\tSCB\tHE\tALB\tHINTS\tPROB_GERMLINE\tPROB_SOMATIC\tPROB_ARTIFACT\tPROB_FFPE_ARTIFACT\tPROB_ABSENT"
  bcftools query -s ffpetumor -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%AF\t%DP\t%SAOBS\t%SROBS\t%OOBS\t%SB\t%ROB\t%RPB\t%SCB\t%HE\t%ALB\t%HINTS]\t%PROB_GERMLINE\t%PROB_SOMATIC\t%PROB_ARTIFACT\t%PROB_FFPE_ARTIFACT\t%PROB_ABSENT\n' \
    ${output_folder}/${id}/data/calls_filt.vcf
) > ${output_folder}/${id}/data/calls_filt.tsv

if [[ "$scenario" == "ffpe" ]]; then
  (
    echo -e "seqnames\tposition\tref\talt\tAF"
    bcftools query -s normal -f '%CHROM\t%POS\t%REF\t%ALT[\t%AF]\n' \
      ${output_folder}/${id}/data/calls_filt.vcf
  ) > ${output_folder}/${id}/data/calls_filt_normal.tsv
fi

python preprocess.py --id $id --vcf $vcf --tumor_cell_quantity $tumor_cell_quantity --output_folder ${output_folder}
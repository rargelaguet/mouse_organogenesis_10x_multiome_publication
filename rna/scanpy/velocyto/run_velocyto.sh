
# velocyto run10x -m repeat_msk.gtf mypath/sample01 somepath/refdata-cellranger-mm10-1.2.0/genes/genes.gtf

indir="/bi/group/reik/ricard/data/gastrulation_multiome_10x/original"

# samples=( "E7.5_rep1" "E7.5_rep2" "E8.0_rep1" "E8.0_rep2" "E8.5_rep1" "E8.5_rep2" "E8.75_rep1" "E8.75_rep2" )
# samples=( "E7.5_rep2" "E8.5_rep1" "E8.5_rep2" "E8.75_rep1" "E8.75_rep2" )
samples=( "E7.75_rep1" "E8.5_CRISPR_T_KO" "E8.5_CRISPR_T_WT" )

threads=1
# mem=1000

mask_file="/bi/group/reik/ricard/data/mm10_sequence/repeats/mm10_rmsk.gtf"

for i in "${samples[@]}"; do
	echo "$i"
	# cmd="velocyto run10x -m ${mask_file} --samtools-threads $threads --samtools-memory 40000 ${indir}/${i} /bi/scratch/Stephen_Clark/annotations/gtf/Mus_musculus.GRCm38.98.gtf"
	cmd="velocyto run10x -m ${mask_file} ${indir}/${i} /bi/scratch/Stephen_Clark/annotations/gtf/Mus_musculus.GRCm38.98.gtf"
	echo $cmd
	sbatch -n $threads --mem 90G --wrap $cmd
done

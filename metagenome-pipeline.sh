##读取宏基因组样本，进行质控、组装、去冗余。
cat ID|while a;do
trimmomatic PE  00-sample/${a}_1.fastq.gz 00-sample/${a}_2.fastq.gz 01-sample/${a}_1.fq.gz unpaired/${a}_1.fq.gz 01-sample/${a}_2.fq.gz unpaired/${a}_2.fq.gz   -phred33 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
bowtie2 -p 8 -x hg38 -1 01-sample/${a}_1.fq.gz -2 01-sample/${a}_2.fq.gz -S 02-bam/${a}.sam
samtools view -bSh 02-bam/${a}.sam > 02-bam/${a}.bam
samtools sort -@ 8  02-bam/${a}.bam -o  02-bam/${a}.sorted.bam
samtools index  02-bam/${a}.sorted.bam
samtools view -b -f 12 -F 256 02-bam/${a}.sorted.bam >  02-bam/${a}.unmapped.bam
samtools sort -n  02-bam/${a}.unmapped.bam  -O BAM -o  02-bam/${a}.unmapped.sort.bam
samtools fastq -@ 8  02-bam/${a}.unmapped.sort.bam -1 03-removeHuman/${a}_1.fq.gz -2 03-removeHuman/${a}_2.fq.gz -n


megahit -1 03-removeHuman/${a}_1.fq -2 03-removeHuman/${a}_2.fq -o 04-megahit/${a}
prokka  04-megahit/${a}/${a}.final.contigs.fa   --outdir 05-prokka/${a} --prefix ${a}  --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses
cd-hit -i 05-prokka/${a}/${a}.ffn -o 06-cd-hit/${a}.nucleotide.fa -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -g 1
transeq -sequence  06-cd-hit/${a}.nucleotide.fa -outseq 06-protein/${a}.protein.fa -trim Y
sed -i 's/_1 / /' 06-protein/${a}.protein.fa
kraken2 --threads 90 --paired --db kraken2_dbtest --report  07-kraken2/${a}.kreport --output 07-kraken2/${a}.kraken 00-sample/${a}_1.fq 00-sample/${a}_2.fq
bracken -d /home/xieyong/Metageome/kraken2_dbtest -i 07-kraken2/${a}.kreport -o 07-kraken2/${a}.bracken.S  -l S -t 90
done

##使用ARGs-OAP管道注释抗生素耐药性基因
args_oap stage_one -i 06-cd-hit -o 06-args  -f fa -t 8
args_oap stage_two -i 08-args -t 8


##使用arg_ranker管道注释抗生素耐药性基因风险水平
arg_ranker -i $INPUT -kkdb kraken2_dbtest
sh arg_ranking/script_output/arg_ranker.sh



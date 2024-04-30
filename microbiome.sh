##读取宏基因组样本，进行质控、组装、去冗余。
cat ID|while a;do
trimmomatic PE  00-sample/${a}_1.fastq.gz 00-sample/${a}_2.fastq.gz 01-sample/${a}_1.fq.gz unpaired/${a}_1.fq.gz 01-sample/${a}_2.fq.gz unpaired/${a}_2.fq.gz   -phred33 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50
megahit -1 01-sample/${a}_1.fq -2 01-sample/${a}_2.fq -o 02-megahit/${a}
prokka  02-megahit/${a}/${a}.final.contigs.fa   --outdir 03-prokka/${a} --prefix ${a}  --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses
cd-hit -i 03-prokka/${a}/${a}.ffn -o 04-cd-hit/${a}.nucleotide.fa -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -g 1
transeq -sequence  04-cd-hit/${a}.nucleotide.fa -outseq 05-protein/${a}.protein.fa -trim Y
sed -i 's/_1 / /' 05-protein/${a}.protein.fa
done

##使用ARGs-OAP管道注释抗生素耐药性基因
args_oap stage_one -i 04-cd-hit -o 06-args  -f fa -t 8
args_oap stage_two -i 06-args -t 8


##含有BacA的contigs和其他样本的比对
#构建数据库
makeblastdb -in  ../06-args/extracted.fa  -dbtype nucl -parse_seqids -out ./index

cat ID|while read a;do
blastx -query ${a}.baca.fa   -db ./index -evalue 1e-5 -outfmt 6 -num_threads 6 -out ${a}-out-file
done

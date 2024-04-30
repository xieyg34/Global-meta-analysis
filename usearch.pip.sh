##安装vsearch,cutadapt,usearch
mkdir test
mkdir temp
mkdir result
mkdir 02-merged
mkdir 03-filt

##合并单双端序列及切除引物与质控

cat V3-PE.id|while read a;do
usearch11 -fastq_mergepairs 01-data/${a}_1.fastq -reverse 01-data/${a}_2.fastq -fastqout test/${a}.merged.fq
cutadapt -g GTGYCAGCMGCCGCGGTAA -o test/${a}-cut1.fq test/${a}.merged.fq
cutadapt -b ATTAGAWACCCVNGTAGTCC -o test/${a}-cut2.fq test/${a}-cut1.fq
usearch10 -fastx_relabel test/${a}-cut2.fq -fastqout 02-merged/${a}.merged.fq -prefix ${a}.
vsearch --fastx_filter 02-merged/${a}.merged.fq  --fastq_maxee_rate 0.01 --fastq_qmax 44  --fastq_maxns 0  --fastaout 03-filt/${a}.filt-fa
done 

cat V3-SE.id|while read a ;do
cutadapt -g GTGYCAGCMGCCGCGGTAA -o test/${a}-cut1.fq 01-data/${a}_1.fastq
cutadapt -b ATTAGAWACCCVNGTAGTCC -o test/${a}-cut2.fq test/${a}-cut1.fq
usearch10 -fastx_relabel test/${a}-cut2.fq -fastqout 02-merged/${a}.merged.fq -prefix ${a}.
vsearch --fastx_filter 02-merged/${a}.merged.fq  --fastq_maxee_rate 0.01 --fastq_qmax 44  --fastq_maxns 0  --fastaout 03-filt/${a}.filt-fa
done

cat V4-PE.id|while read a;do
cutadapt -g GTGYCAGCMGCCGCGGTAA -G GGACTACNVGGGTWTCTAAT -o test/${a}_cut1.fq -p test/${a}_cut2.fq  01-data/${a}_1.fastq  01-data/${a}_2.fastq
vsearch --fastq_mergepairs test/${a}_cut1.fq --reverse test/${a}_cut2.fq   --fastqout 02-merged/${a}.merged.fq --relabel ${a}. 
vsearch --fastx_filter 02-merged/${a}.merged.fq  --fastq_maxee_rate 0.01 --fastq_qmax 44  --fastq_maxns 0  --fastaout 03-filt/${a}.filt-fa
done


cat V4-SE.id|while read a;do
cutadapt -g GTGYCAGCMGCCGCGGTAA -o test/${a}-cut1.fq 01-data/${a}_1.fastq
cutadapt -b ATTAGAWACCCVNGTAGTCC -o test/${a}-cut2.fq test/${a}-cut1.fq
usearch10 -fastx_relabel test/${a}-cut2.fq -fastqout 02-merged/${a}.merged.fq -prefix ${a}.
vsearch --fastx_filter 02-merged/${a}.merged.fq  --fastq_maxee_rate 0.01 --fastq_qmax 44  --fastq_maxns 0  --fastaout 03-filt/${a}.filt-fa
done

## 合并所有样品至同一文件

cat 03-filt/*.filt-fa > temp/filtered.fa


## 序列去冗余
time vsearch --derep_fulllength temp/filtered.fa  --output temp/uniques.fa --relabel Uni --minuniquesize 10 --sizeout
##生成OTU 97%聚类OTU
time usearch10 -cluster_otus temp/uniques.fa -otus temp/otus.fa -relabel OTU_

## vsearch+silva去嵌合
time vsearch --uchime_ref temp/otus.fa   --db ../db/silva_16s_v123.fa  --nonchimeras result/otus.fa
## 生成OTU表
time vsearch --usearch_global temp/filtered.fa --db result/otus.fa   --otutabout result/otutab.txt --id 0.97 --threads 4

## OTU表简单统计 Summary OTUs table
usearch10 -otutab_stats result/otutab.txt  -output result/otutab.stat
## 等量抽样标准化
usearch10 -otutab_norm result/otutab.txt   -sample_size 30000  -output result/otutab_norm.txt
##物种注释
usearch10 -sintax result/otus.fa -db db/rdp_16s_v16_sp.fa  -strand both -tabbedout result/sintax.txt -sintax_cutoff 0.6
cut -f 1,4 result/sintax.txt  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//'   > result/taxonomy2.txt

# 生成物种表格：注意OTU中会有末知为空白，补齐分类未知新物种为Unassigned
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Una
ssigned";a["
s"]="Unassigned"; split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];}   print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}'   
result/taxonomy2.txt >temp/otus.tax

sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/'  > result/taxonom
y.txt
head -n3 result/taxonomy.txt

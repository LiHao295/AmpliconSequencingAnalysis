##########################################
#####      Summarized by Li Hao      #####
#####           20180117             #####
#####      lihao@stu.xmu.edu.cn      #####
#####        usearch，mothur         #####
##### http://www.drive5.com/usearch/ #####
#####    https://www.mothur.org/     #####
##########################################

###1.Make contigs, filter the seqs and split groups by mothur
mothur "#make.contigs(ffastq=L1_1.clean.fq,rfastq=L1_2.clean.fq,processors=50);trim.seqs(fasta=L1_1.clean.trim.contigs.fasta,oligos=V4F.txt,minlength=200,maxlength=275,maxambig=0,maxhomop=8,bdiffs=0,pdiffs=0,tdiffs=0);trim.seqs(fasta=L1_1.clean.trim.contigs.scrap.fasta,oligos=V4R.txt,minlength=200,maxlength=275,maxambig=0,maxhomop=8,bdiffs=0,pdiffs=0,tdiffs=0,flip=T);"
cat hos_L1_1.clean.trim.contigs.scrap.trim.fasta hos_L1_1.clean.trim.contigs.trim.fasta > hos.fasta
cat hos_L1_1.clean.trim.contigs.groups hos_L1_1.clean.trim.contigs.scrap.groups > hos.group
mothur "#split.groups(fasta=hos.fasta,group=hos.group,groups=all);"
 
###2.Data preperation
 rename 's/hos./hos_/' *
 fa=`ls pwd|grep .fasta$`
 for i in $fa;do k=`basename $i .fasta`;perl rename_fasta.pl -i $i -o $i.fa -n $k;done
 cat *.fa >hos_v4.fasta

###3.Analysis by USEARCH
###3.1  unique the seqs
usearch -fastx_uniques hos_v4.fasta -fastaout hos_v4_unique.fasta -minuniquesize 1 -sizeout -threads 30

###3.2  remove the chimera and make Zotu or make 97% OTUs and filter chimeras 
usearch -unoise3 hos_v4_unique.fasta  -minsize 8 -zotus zotus.fa
#usearch -cluster_otus uniques.fa -otus otus.fa -relabel Otu

###3.4 Make OTU table 
sed 's/>Zotu/>Otu/g' zotus.fa >zotu.fa
usearch -otutab hos_v4.fasta -zotus zotu.fa -otutabout zotutab_raw.txt -mothur_shared_out zotutab_raw4mothur.txt -threads 55

###3.5 remove the seqs by coss talking
usearch -uncross zotutab_raw.txt -tabbedout out.txt -report rep.txt -otutabout zotutab_nocosstalk.txt

###3.6 Normalize to 5k reads / sample
usearch -otutab_norm zotutab_nocosstalk.txt -sample_size 5000 -output zotutab.txt

###3.7 Alpha diversity
usearch -alpha_div otutab.txt -output alpha_div

###3.8 Make OTU tree
usearch -cluster_agg zotu.fa -treeout otus.tree

###3.9 Beta diversity
mkdir beta/
usearch -beta_div zotutab.txt -tree otus.tree -filename_prefix beta/

###
Rscript /data/LH/amplicon/usearch/trans_loewr_triangle.R ./beta/unifrac.sorted.txt unifrac.sorted2.txt
mothur "#nmds(phylip=unifrac.sorted_lt.txt);pcoa(phylip=current)"

###3.10 Predict taxonomy by mothur
mothur1.36 "#classify.seqs(fasta=zotu.fa,reference=eztaxon_qiime_full.fasta, taxonomy=eztaxon_id_taxonomy.tax, cutoff=80,processors=50)"



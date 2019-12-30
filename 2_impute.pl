#£¡/usr/bin/perl -w
use strict;

## remove low quality SNPs and samples
system("plink --noweb --bfile ../AMP-AD_HBTRC_MSSM_IlluminaHumanHap650Y --exclude remove_snp.txt --remove remove_patients.txt  --make-bed --out input2");

## split the genotyping data into 22 chromomsome
for(my $pp=1; $pp < 23; $pp++){
  system("plink --noweb --bfile input2 --chr $pp --make-bed --out data/input_chr$pp"); 
}

## estimation of haplotypes
for(my $pp=1; $pp < 23; $pp++){
  system("shapeit --force --input-bed data/input_chr$pp.bed data/input_chr$pp.bim data/input_chr$pp.fam --input-map ~/work/impute/1000G/genetic_map_chr${pp}_combined_b37.txt --thread 45 --effective-size 11418 --output-max data/chr$pp.haps data/chr$pp.sample --output-log data/chr$pp.sample.log");
}


## impute
for(my $pp=1; $pp < 23; $pp++){
  system("gunzip -cd ~/work/impute/1000G/1000GP_Phase3_chr$pp.legend.gz > 1000.legend");
  my $xx=`tail -n 1 1000.legend`;
  my @xx=split(/\s/,$xx);
  my $max=$xx[1];
  my $zz=0;
  my $step=3000000;
  for(my $i=250000;$i< $max;$i=$i+$step){
    open OUT,">jobs/job${pp}_$zz.sh" or die;
    my $j=$i+$step;
    print OUT "impute2 -filt_rules_l 'EUR<0.01' -m ~/work/impute/1000G/genetic_map_chr${pp}_combined_b37.txt -known_haps_g data/chr$pp.haps -h ~/work/impute/1000G/1000GP_Phase3_chr$pp.hap.gz  -l ~/work/impute/1000G/1000GP_Phase3_chr$pp.legend.gz -int $i $j  -Ne 20000  -o result/imputed_chr${pp}_$zz.gen -seed 367946\n";
    $zz++;
    close OUT;
  } 
}

## run them on a cluster system
my @files=<jobs/*.sh>;    
foreach my $file (@files){
    system("qsub $file &");
}
## GEN to PED conversion
for(my $pp=1; $pp < 23; $pp++){
  system("cat result/imputed_chr${pp}_*.gen > result/imputed_chr${pp}.gen");
  system("nohup gtool -G --g result/imputed_chr${pp}.gen --s data/chr$pp.sample --ped result/impute2_chr${pp}.ped --map result/impute2_chr${pp}.map --threshold 0.9 &");
}



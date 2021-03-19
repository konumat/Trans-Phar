

###echo inputted ICD10 code
echo "ICD10 code" $1


#####munge sumstats for TWAS
cd ../Input_GWASsummary

for filename in $( ls *.sumstats ); do

echo ${filename}

mkdir -p ../Output/${filename::-9}

../focus/bin/focus munge ${filename} --output ../Output/${filename::-9}/${filename::-9}_munge


#####TWAS

cd ../Output/${filename::-9}
mungedfile=$( ls *.gz )
mkdir -p ./TWASresults

  cd ../../db
  for tissuename in $( ls *.db ); do

  mkdir -p ../Output/${filename::-9}/TWASresults/${tissuename::-3}

    for i in {1..22}; do ../focus/bin/focus finemap ../Output/${filename::-9}/${mungedfile} ../1000G_EUR_Phase3_plink/1000G.EUR.QC.$i ${tissuename} --chr $i --tissue ${tissuename:5:-3} --p-threshold 1  --out ../Output/${filename::-9}/TWASresults/${tissuename::-3}/${tissuename::-3}.chr$i; done
    cd ../Output/${filename::-9}/TWASresults/${tissuename::-3}
    Rscript ../../../../script/FOCUS_assoc_test_unifyingresult.R
    python3 ../../../../script/FOCUS_resultstoplots.py
    ##特に必要なファイルは別フォルダに
    mkdir -p ../ALLTISSUE
    mv *.png ../ALLTISSUE
    mv *chr_all.focus.tsv ../ALLTISSUE
    cd ../../../../db


  done

#####Spearman's negative correlation analysis (all cell-type - drug pairs in each tissue/cell-type category)

cd ../Output/${filename::-9}
mkdir -p ./Spearmanresults

Rscript ../../script/Spearman.R $1


#####After ALL procedures, input GWAS summary statistics were transferred to other folder

mv ../../Input_GWASsummary/${filename} ../../Input_GWASsummary_done/${filename}

cd ../../Input_GWASsummary

done
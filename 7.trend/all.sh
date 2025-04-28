#Trend pipe start now!
#compute-0-57.local
mkdir -p /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/
cp /Bio/Bin/pipeline/Small_pipe/v4.5/trend/config/Trend /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
mkdir -p /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/
cp /Bio/Bin/pipeline/Small_pipe/v4.5/trend/config/readme.docx /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/
mkdir -p /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/result
mkdir -p /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/enrich
#Format input exp data : instead of [0,-]
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/format_exp.pl /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/exp.xls l /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/rpkm.xls
#Run stem
sed -i 's/Data_File	aaa/Data_File	\/kfs\/Bio\/Project\/202303\/GHI22120213-1\/GHI22120213-1_sup_11\/trend\/out\/rpkm.xls/' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
sed -i 's/K-means]	aaa/K-means]	STEM Clustering Method/' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
sed -i 's/Minimum_Absolute_Log_Ratio_Expression	aaa/Minimum_Absolute_Log_Ratio_Expression	1.0/' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
sed -i 's/Maximum_Number_of_Model_Profiles	aaa/Maximum_Number_of_Model_Profiles	20/' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
sed -i 's/Maximum_Unit_Change_in_Model_Profiles_between_Time_Points	aaa/Maximum_Unit_Change_in_Model_Profiles_between_Time_Points	1/' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
sed -i 's/add 0]	aaa/add 0]	Log normalize data/' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput/Trend
export DISPLAY=:1 && /Bio/Bin/pipeline//System_Programs/java -mx1024M -jar /Bio/Bin/pipeline//System_Programs/my_SoftWare/trend/stem.jar -b /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/steminput  /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/
#Draw trend
sed -i 's/,//g' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.txt
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/format_trend_profiletable.pl /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_profiletable.txt /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_profiletable.txt2 0
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/format_genetable.pl /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.txt /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.txt2 0
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/R_plot.pl Trend_profile_table -trend_genetable /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.txt2 -trend_profiletable /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_profiletable.txt2 -method log -pvalue 0.05 -name_key Genes -profile_name profile -outdir /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/R_plot.pl Trend_profile_fig -pvalue 0.05 -infile /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_profiletable.txt2 -name_key Genes -profile_name profile -outdir /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/R_plot.pl Trend_profile_line -name_key Genes -profile_name profile -infile /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_profiletable.txt2 -trend_genetable /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.txt2 -outpfx /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/Trend_profile_line
#finish draw profile picture
#Split profile
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/changeid.pl /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/exp.xls /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.txt2 >/kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.xls
cp /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.xls /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/all_profile.xls
#Add annot and get glist
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/Get_Profile_Glist.pl -code utf8 -profile_name profile /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_genetable.xls /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/exp.xls 
#Draw Profile
sed 's/Gene/Gene/g' /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/stemoutput/Trend_profiletable.txt2 >/kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/profile_stat.xls
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/R_plot.pl Trend_rect_profile -infile /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/profile_stat.xls  -outpfx /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/profile_stat -title Genes_in_Profile -header yes -profile_name profile
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/enrich.pl -i /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/enrich -ref /public/Database/pipRef/Homo_sapiens/Ensembl_release93/ref4SingleCell/GRCh38_annot/annot/ref -run_ko 1 -run_go 1 -run_do 1 -run_reactome 1 -type nodiff -key gene
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/path_sta.pl -i "/kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/enrich/KO/*.path.xls" -o /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/enrich/KO/all_pathway.xls -map_title /public/Database/kegg_hsm/20220220/map_title.tab
#make link
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/make_link.pl /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/Trend_analysis/ /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/result yes
#write report
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/report/trend_report.pl -indir /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/result 1 1 -enrich 1 -ratio 1.0 -exp /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/exp.xls -run_go 1 -run_do 1 -run_reactome 1 -met 0 -nonlazy 1
perl /Bio/Bin/pipeline//System_Programs/Mono_Report_Generator/v2.1/Mono_Report_Generator.pl /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/result/index.html >/kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/GHI22120213-1_sup_8—Epi.html
perl /Bio/Bin/pipeline/Small_pipe/v4.5/trend/src/report/trend_report.pl -indir /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out/result 1 1 -enrich 1 -ratio 1.0 -exp /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/exp.xls -run_go 1 -run_do 1 -run_reactome 1 -met 0 -nonlazy 0
#tar result
cd /kfs/Bio/Project/202303/GHI22120213-1/GHI22120213-1_sup_11/trend/out;ln -sf --no-dereference result GHI22120213-1_sup_8—Epi;tar jcf GHI22120213-1_sup_8—Epi.tar.bz2 --dereference GHI22120213-1_sup_8—Epi
#All Done~~~~

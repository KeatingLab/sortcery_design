
for each in `cat dirs`;
	do 
	cd $each/workspace/;
	for bar in barcode*;
		do 
		sed -ne '1~8p;2~8p;3~8p;4~8p' $bar > /scratch/users/vxue/data/vxue_mit_edu/$each'_'$bar'_'0;
		sed -ne '5~8p;6~8p;7~8p;8~8p' $bar > /scratch/users/vxue/data/vxue_mit_edu/$each'_'$bar'_'1
	done;
	cd ../../;
done;


#sed -ne '1~8p;2~8p;3~8p;4~8p' x.fastq > x_1.fastq
#sed -ne '5~8p;6~8p;7~8p;8~8p' x.fastq > x_2.fastq

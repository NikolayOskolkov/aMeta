conda activate gargammel
cd /home/nikolay/WABI/A_Gotherstrom/gargammel/gargammel

endoFrac=(0.1 0.15 0.2 0.3 0.35 0.4 0.5 0.55 0.6 0.65)
contFrac=(0.2 0.15 0.1 0.2 0.15 0.1 0.1 0.15 0.1 0.05)
bactFrac=(0.7 0.7 0.7 0.5 0.5 0.5 0.4 0.3 0.3 0.3)
for i in {0..9}
do
echo \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\# SIMULATING SAMPLE $i WITH GARGAMMEL \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
cut -f1,$[$i+2] ground_truth_ancient_microbes.txt > data_ancient/bact/list
cut -f1,$[$i+2] ground_truth_modern_microbes.txt > data_modern/bact/list
gargammel -n 250000 --comp ${bactFrac[$i]},${contFrac[$i]},${endoFrac[$i]} -rl 125 --loc 3.7424069808 --scale 0.2795148843 -damage 0.03,0.4,0.01,0.3 -o data_ancient/simulation data_ancient/
gargammel -n 250000 --comp ${bactFrac[$i]},${contFrac[$i]},${endoFrac[$i]} -rl 125 -f src/sizefreq.size.gz -o data_modern/simulation data_modern/
cat data_ancient/simulation_s1.fq.gz data_modern/simulation_s1.fq.gz > data/sample$[${i}+1].fastq.gz
done


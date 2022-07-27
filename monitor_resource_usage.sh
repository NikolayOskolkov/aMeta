i=0
rm cpu_usage.txt
rm memory_usage.txt
while true
do
ps -u nikolay -eo %mem,%cpu --sort=-%mem | awk '{ hr=$1 ; sum +=hr} END {print sum}' >> memory_usage.txt
ps -u nikolay -eo %mem,%cpu --sort=-%cpu | awk '{ hr=$2 ; sum +=hr} END {print sum}' >> cpu_usage.txt
sleep 60
i=$[$i+1]
echo The process has been running for ${i} min
done

df<-as.numeric(readLines("memory_usage.txt"))
df<-(df/100)*512
plot(df,type="l",xlab="TIME (MIN)",ylab="MEMORY (RAM) USAGE (GB)")

df_cpu<-as.numeric(readLines("cpu_usage.txt"))
plot(df_cpu,type="l",xlab="TIME (MIN)",ylab="CPU USAGE (%)")



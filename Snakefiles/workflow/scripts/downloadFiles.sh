#! bin/bash 

while read p;
do
	wget ${p} -P $2
done < $1

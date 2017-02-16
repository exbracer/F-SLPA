# !/bin/bash

# input & and output file location
DATASET_LOCATION=../input/
RESULT_LOCATION=../result/evaluation/

# parameters
MAX_T=1000

# for magellan
DATASET=(
	com-amazon.data
	com-dblp.data
	com-lj.data
)

NUM_THREADS=(
	1
	2
	4
	8
	12
	16
	24
	32
	48
)

VERSION=(
	11
	13
	15
	17
	22
	23
	24
	25
	26
)

# executable file
EXECUTABLE_FILE=(
	SLPA
	SLPA_tc
)

# const
UNDERLINE=_
datasetPath=
maxT=$MAX_T
LOG_FILE=runtime.log
echo START!

for exec in ${EXECUTABLE_FILE[*]}
do
	for dataset in ${DATASET[*]}
	do
		datasetPath=$DATASET_LOCATION
		datasetPath+=$dataset
		for version in ${VERSION[*]}
		do
			for numThreads in ${NUM_THREADS[*]}
			do
				resultPath=$RESULT_LOCATION$exec$UNDERLINE$dataset$UNDERLINE$version$UNDERLINE$numThreads

				./$exec -i $datasetPath -t $maxT -M $numThreads -v $version > $LOG_FILE

				cat $LOG_FILE
				resultPath1=$resultPath
				resultPath1+=$UNDERLINE
				resultPath1+=RUNTIME_LOG

				mv ./$LOG_FILE $resultPath1
			done
		done
		./$exec -i $datasetPath -t $maxT -M 1 -v 10 > $LOG_FILE
		cat $LOG_FILE
		resultPath1=$RESULT_LOCATION$exec$UNDERLINE$dataset$UNDERLINE
		resultPath1+=10
		resultPath1+=$UNDERLINE
		resultPath1+=1
		resultPath1+=$UNDERLINE
		resultPath1+=$LOG_FILE
		mv ./$LOG_FILE $resultPath1

		./$exec -i $datasetPath -t $maxT -M 1 -v 20 > $LOG_FILE
		cat $LOG_FILE
		resultPath1=$RESULT_LOCATION$exec$UNDERLINE$dataset$UNDERLINE
		resultPath1+=20
		resultPath1+=$UNDERLINE
		resultPath1+=1
		resultPath1+=$UNDERLINE
		resultPath1+=$LOG_FILE
		mv ./$LOG_FILE $resultPath1
	done
done
#./SLPA -i ../input/com-amazon.data -t 1000 -M 1 -v 17
#

















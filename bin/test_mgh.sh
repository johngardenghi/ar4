#!/bin/bash

param=""

for i in $*
do
    param="$param $i"
done

echo $param

for i in {1..35}
do
    ./mgh $param -prob $i -pl 0
    # ./mgh -p 2 -prob $i
done

if [ -e "results_tab" ]
then
    echo "" >> results_tab
    echo "" >> results_tab
    echo "" >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo "" >> results_tab
    echo "" >> results_tab
fi

if [ -e "results_tab_p2" ]
then
    echo "" >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo "                                                                            p=2                                                                            " >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo "                                         PROBLEM    n    m | FLAG     OBJ FUNCTION     ||GRADIENT||    #ITER   #EVALF   #EVALG       TIME      SIGMA NSTEP " >> results_tab
    cat results_tab_p2 >> results_tab
fi

if [ -e "results_tab_p3" ]
then
    echo "" >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo "                                                                            p=3                                                                            " >> results_tab
    echo " ==========================================================================================================================================================" >> results_tab
    echo "                                         PROBLEM    n    m | FLAG     OBJ FUNCTION     ||GRADIENT||    #ITER   #EVALF   #EVALG       TIME      SIGMA NSTEP " >> results_tab
    cat results_tab_p3 >> results_tab
fi

rm -Rf results_tab_*

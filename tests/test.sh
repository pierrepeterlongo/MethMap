#!/bin/bash
python ../MethMap.py data_test/bank.fa data_test/query.fa 12 False > current_non_converted.afac
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on command python ../MethMap.py data_test/bank.fa data_test/query.fa 12 False"
rm -f *afac*
exit 1
fi

diff current_non_converted.afac res_non_converted.txt
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on diff non converted"
rm -f *afac*
exit 1
fi


python ../MethMap.py data_test/bank.fa data_test/query.fa 12 True > current_converted.afac
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on command python ../MethMap.py data_test/bank.fa data_test/query.fa 12 True"
rm -f *afac*
exit 1
fi

diff current_converted.afac res_converted.txt
if [ $? -ne 0 ] ; then
echo "*** Test: FAILURE on diff converted"
rm -f *afac*
exit 1
fi


echo "*** Test: OK"
rm -f *afac*
exit 0
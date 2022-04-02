#!/usr/bin/env python
# coding: utf-8

# # splits XMFA file for separate input for clonalframe  


import os

in_xmfa = '/ebio/ag-mccann/projects/rso-mar/code/geneTrack/mar1_concat.xmfa'
#'/Users/devseeva/Desktop/work/sm_workflow/geneTrack/phy2_test_outputs/phy2_test_concat.xmfa'

l_num = int(os.popen("grep -c = "+in_xmfa).read())
split_num = round(l_num / 2)
split_num


xmfa = open(in_xmfa, 'r')
out1 = open(in_xmfa.replace('.xmfa', '_s1.xmfa'), 'w')
out2 = open(in_xmfa.replace('.xmfa', '_s2.xmfa'), 'w')

counter = 0
for l in xmfa:
    if counter < split_num:
        out1.write(l)
    else:
        out2.write(l)
    if '=' in l:
        counter = counter + 1

print(counter)    
xmfa.close()
out1.close()
out2.close()

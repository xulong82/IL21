#!/bin/bash

cd $HOME/Dropbox/GitHub/Lupus/homer

$HOME/Applications/homer/bin/findMotifs.pl N.txt mouse N/ -p 6 
$HOME/Applications/homer/bin/findMotifs.pl ACT.txt mouse ACT/ -p 6 
$HOME/Applications/homer/bin/findMotifs.pl IL21.txt mouse IL21/ -p 6 


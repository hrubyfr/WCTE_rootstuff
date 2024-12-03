#!/bin/bash
rsync -avzP -e "ssh -p 2233" mpmt@192.168.11.252:~/ReadoutTest_Oct2024/WCTEDataFormatExample/tuples ~/Desktop

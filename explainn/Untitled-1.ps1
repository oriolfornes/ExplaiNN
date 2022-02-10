$ zless ./results/IRF4/CAM/ENCODE/IRF4-BATF/IRF4.WT.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$aice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$aice+=1;}}print("AICE = $aice\nIRF4 = $irf4\n");'
AICE = 12980
IRF4 = 15277
$ zless ./results/IRF4/CAM/ENCODE/IRF4-BATF/IRF4.T95R.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$aice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$aice+=1;}}print("AICE = $aice\nIRF4 = $irf4\n");'
AICE = 27790
IRF4 = 14993
$ zless ./results/IRF4/CAM/ENCODE/IRF4-BATF/IRF4.WT-unique.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$aice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$aice+=1;}}print("AICE = $aice\nIRF4 = $irf4\n");'
AICE = 5210
IRF4 = 10759
$ zless ./results/IRF4/CAM/ENCODE/IRF4-BATF/IRF4.T95R-unique.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$aice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$aice+=1;}}print("AICE = $aice\nIRF4 = $irf4\n");'
AICE = 19765
IRF4 = 10732
$ zless ./results/IRF4/CAM/ENCODE/IRF4-SPI1/IRF4.WT.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$eice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$eice+=1;}}print("EICE = $eice\nIRF4 = $irf4\n");'
EICE = 8787
IRF4 = 19470
$ zless ./results/IRF4/CAM/ENCODE/IRF4-SPI1/IRF4.T95R.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$eice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$eice+=1;}}print("EICE = $eice\nIRF4 = $irf4\n");'
EICE = 8601
IRF4 = 34182
$ zless ./results/IRF4/CAM/ENCODE/IRF4-SPI1/IRF4.T95R-unique.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$eice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$eice+=1;}}print("EICE = $eice\nIRF4 = $irf4\n");'
EICE = 4783
IRF4 = 25714
$ zless ./results/IRF4/CAM/ENCODE/IRF4-SPI1/IRF4.WT-unique.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$eice=0;$irf4=0;while(<>){chomp;if($_<.5){$irf4+=1;}else{$eice+=1;}}print("EICE = $eice\nIRF4 = $irf4\n");'
EICE = 4925
IRF4 = 11044
(cam) oriol@GPURTX:/mnt/md1/home/oriol/CAM$ zless ./results/IRF4/CAM/ChIP-seq/WT-T95R-unique/IRF4-BATF.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$t95r=0;$wt=0;while(<>){chomp;if($_<.5){$wt+=1;}else{$t95r+=1;}}print("T95R = $t95r\nWT = $wt\n");'
T95R = 6346
WT = 7541
(cam) oriol@GPURTX:/mnt/md1/home/oriol/CAM$ zless ./results/IRF4/CAM/ChIP-seq/WT-T95R-unique/IRF4-noBATF.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$t95r=0;$wt=0;while(<>){chomp;if($_<.5){$wt+=1;}else{$t95r+=1;}}print("T95R = $t95r\nWT = $wt\n");'
T95R = 2009
WT = 6379
(cam) oriol@GPURTX:/mnt/md1/home/oriol/CAM$ zless ./results/IRF4/CAM/ChIP-seq/WT-T95R-unique/IRF4-SPI1.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$t95r=0;$wt=0;while(<>){chomp;if($_<.5){$wt+=1;}else{$t95r+=1;}}print("T95R = $t95r\nWT = $wt\n");'
T95R = 2219
WT = 6252
(cam) oriol@GPURTX:/mnt/md1/home/oriol/CAM$ zless ./results/IRF4/CAM/ChIP-seq/WT-T95R-unique/IRF4-noSPI1.tsv.gz | grep -v SeqId | cut -f 4 | perl -e '$t95r=0;$wt=0;while(<>){chomp;if($_<.5){$wt+=1;}else{$t95r+=1;}}print("T95R = $t95r\nWT = $wt\n");'
T95R = 5885
WT = 7478
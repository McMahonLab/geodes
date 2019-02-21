#!/bin/bash
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,$30,$31,$32,$33,$34,$35}' GEODES_ID90_2018-03-01.txt > Sparkling_ID90_2018-03-02.txt
head -1 Sparkling_ID90_2018-03-02.txt > colnames.txt
awk 'NR > 1{s=0; for (i=3;i<=NF;i++) s+=$i; if (s!=0)print}' Sparkling_ID90_2018-03-02.txt > temp.txt
cat colnames.txt temp.txt > Sparkling_ID90_2018-03-02.readcounts.txt


# Trout Bog (epilimnion) - GEODEs053-100, 110
awk '{print $1,$36,$37,$38,$39,$40,$41,$42,$43,$44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56, $57}' GEODES_ID90_2018-03-01.txt > Trout_ID90_2018-03-02.txt
head -1 Trout_ID90_2018-03-02.txt > colnames.txt
awk 'NR > 1{s=0; for (i=3;i<=NF;i++) s+=$i; if (s!=0)print}' Trout_ID90_2018-03-02.txt > temp.txt
cat colnames.txt temp.txt > Trout_ID90_2018-03-02.readcounts.txt



# Mendota - GEODES113-164
awk '{print $1,$58,$59,$60,$61,$62,$63,$64,$65,$66,$67,$68,$69,$70,$71,$72,$73,$74,$75,$76,$77,$78,$79,$80,$81,$82,$83,$84,$85,$86,$87,$88}' GEODES_ID90_2018-03-01.txt > Mendota_ID90_2018-03-02.txt
head -1 Mendota_ID90_2018-03-02.txt > colnames.txt
awk 'NR > 1{s=0; for (i=3;i<=NF;i++) s+=$i; if (s!=0)print}' Mendota_ID90_2018-03-02.txt > temp.txt
cat colnames.txt temp.txt > Mendota_ID90_2018-03-02.readcounts.txt


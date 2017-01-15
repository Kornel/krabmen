#!/usr/bin/env sh

awk -f filter.awk < BRCA.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt > out.txt 

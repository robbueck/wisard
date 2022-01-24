#!/bin/sh
PHAGES=()
for f in /data/meyer/rob/phages/robert_phage/genome_alignments/method_comparison/other_data/HIV/*.fasta;
do PHAGES+=( "$f" );
done
dir=/data/meyer/rob/phages/robert_phage/genome_alignments/method_comparison/other_data/HIV/
for s in ${PHAGES[@]}; do PHAGES=( "${PHAGES[@]/$s}" ) ;
    for q in ${PHAGES[@]}; 
        do n=${q%.fasta};
        n=$(basename "$n")
        m=${s%.fasta};
        m=$(basename "$m")
        o="${m}_${n}" ;
        ~/Programs/ab-blast/ab-blast-20200317-linux-x64/ab-blastn ${s%.fasta} $q \
            M=1 N=-1 Q=3 R=2 W=9 -hspmax 0 -kap -wordmask=seg \
            -mformat=7,"${dir}/${o}_ab_blastn.xml";
        /data/meyer/rob/Programs/ab-blast/ab-blast-20200317-linux-x64/ab-tblastx ${s%.fasta} $q \
            -hspmax=0 -matrix=BLOSUM62 -W=3 -T=50 -hitdist=40  -kap -filter=seg -dbgcode=11 -mformat=7,"${dir}/${o}_ab_tblastx.xml"
        echo $o;
    done;
done

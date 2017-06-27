# bam2circos
Visualize your BAM alignments with Circos

# Genbank identifiers
Circos does not support the vertical bar `|` character in the reference name. Therefore, if your BAM file contains reference names of the format `gi|9630643|ref|NC_001918.1|`, you will have to remove the vertical bars from your BAM file. An easy way to do this is with the following one liner given that your BAM file is called foo.bam and you have samtools version < 1.

```bash
samtools view -h foo.bam | perl -pe 's/gi\|\d+\|ref\|(\S+)\|/$1/' | samtools view -bS - | samtools sort - foo.reformat
```

If you have samtools version >1, change the `samtools sort` command to the following `samtools sort -o foo.reformat.bam`.


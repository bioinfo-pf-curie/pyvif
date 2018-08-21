# PyVIF Usage

## Standalone

For data fast analysis, the standalone generates a standalone HTML report that contains
quality control and list of detected breakpoints. Subreads alignments with `minimap2` on
human and virus genomes are needed for PyVIF.

```
minimap2 -t 8 \
         -La -x map-pb hg38.fasta subreads.fastq \
         | samtools sort -@ 7 \
         -o human.bam \
         && samtools index human.bam

minimap2 -t 8 \
         -La -x map-pb virus.fasta subreads.fastq \
         | samtools sort -@ 7 \
         -o virus.bam \
         && samtools index virus.bam
```

Then, run `pyvif` with those bam files.

```
pyvif -u human.bam -v virus.bam -o pyvif_report.html
```

## API

The second approach is to use the package directly in ipython.
Run the notebook with jupyter `jupyter notebook` and open the notebook `notebooks/pyvif.ipynb`.

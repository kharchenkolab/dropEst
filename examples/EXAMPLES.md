# Pipeline examples

## Quick start: InDrop v3 mouse data
To run the pipeline on the small subset of inDrop Mouse BMCs, use the following commands.

Go to the script directory:
```bash
cd dropEst/examples/scg71_demo
``` 

Download data and prepare structure of the folders:
```bash
make
``` 

Run the pipeline, using STAR aligner. It requires you to have `STAR` and `zcat` in the PATH:
```bash
./pipeline.sh path_to_dropest_binaries path_to_config path_to_star_index path_to_genes_gtf
```

Here:
- `path_to_dropest_binaries` - path to the folder with `dropest` and `droptag` binaries;
- `path_to_config`: path to indrop_v3.xml config;
- `path_to_star_index`: path to STAR index;
- `path_to_genes_gtf`: path to gtf file with annotated genes.

Example:
```bash
./pipeline.sh ~/dropEst/build ~/dropEst/configs/indrop_v3.xml /n/groups/pklab/genomes/mm10/STARIndex /n/groups/pklab/genomes/mm10/genes.gtf
```

The following command downloads `mm10_genes.gtf`, which was used by the authors:
```bash
make genes
```

To examine validity of the results, you can download whole snapshot of the folder, which we got by running the same commands.
The following command downloads them to `results` folder:
```bash
make results
```
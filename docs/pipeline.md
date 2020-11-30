# Pipeline walkthrough


### 1 Organise raw calls into a data directory (here `data-dir`) with the following structure:

```
data-dir/
    lumpy/[.lumpy.gt_all.vcf]
    delly/[.delly.vcf]
    novobreak/[.novoBreak.pass.flt.vcf]
    cnv/
    freec/[_sig_cnvs.txt]
```    

### 2 Organise output from CNV-Seq into a directory (here `cnv-data`) that contains all `.cnv` files for small window run (e.g. 500 bps)


### 3 Run `freecFilt.sh` from `freec` directory to output filtered Control-freec calls:

```bash
cd freec/
bash freecFilt.sh *_sig_cnvs.txt
```
 
 
### 4 Run `run_parser.sh` ([svParser](https://github.com/nriddiford/svParser) wrapper):

```bash
data=/path/to/data-dir
cnv_dir=/path/to/cnv-dir/w500
output_directory=/path/to/write/filtered_data

# -f [filter]
# -m [merge filtered calls]
# -a [annontate filtered calls]
# -s [somatic tumour only]

bash runParser.sh
    -d "${data}" \
    -c "${cnv_dir}" \
    -o "${output_directory}" \
    -fma \
    -s
```
This will generate the following files in the specified output_directory `${output_directory}`:

* Filtered vcfs to `${output_directory}`
* Summary files for each caller (delly/lumpy/novobreak) `${output_directory}/summary/[.filtered.summary.txt]`
* CNV-annotated files for each caller `${output_directory}/summary/[.filtered.summary.cnv.txt]`
* Merged & clustered (per sample) temp files
* Annotated files merged per sample `${output_directory}/summary/merged/[sample_annotated_SVs.txt]`


### 5 Inspect calls in IGV

```bash
script_bin=/path/to/svParser/scripts
small_cnv_calls=${cnv_dir}
large_cnv_calls=/path/to/cnv-dir/w1000

sample_name=sample1
```


#### 5.1 Plotting SVs using [cnvPlotteR](https://github.com/nriddiford/cnvPlotteR):

```python
python "${script_bin}/svPlotter.py" \
    -f "${sample_name}_annotated_SVs.txt" \
    --cnv_small "${small_cnv_calls}/${sample_name}_*.cnv" \
    --cnv_large "${large_cnv_calls}/${sample_name}_*.cnv" \
    -t
```

Copy and paste a variant from the output and use as R command in cnvPlotteR


#### 5.2 Manually adding/adjusting variants:

CNV events can be added to `${sample_name}_annotated_SVs.txt` by providing the following information:


#### 6 Annotate variants with allele frequency (& filter etc) using [svSupport](https://nriddiford.github.io/svSupport)

* Run `bash scripts/run_all.sh -cs -v /path/to/${sample_name}_annotated_SVs.txt/files` to make configs for all variant files in specified path and run svSupport


#### 7 Stitch together complex variants

* Copy output from `svSupport` (`${sample_name}_svSupport.txt`) to svParser output dir `${output_directory}/summary/merged`
```{bash}
for f in *._svSupport.txt
do
  python ${script_bin}/svStitch -w 5500 -
done
```
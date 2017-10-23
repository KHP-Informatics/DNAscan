# DNA-NGS-scan

DNAscan is a fast and efficient bioinformatics pipeline that allows for the analysis of DNA Next Generation sequencing data, requiring very little computational effort and memory usage. DNAscan can analyse raw whole genome NGS data in ~8 hours, using as little as 4 cpus and ~4 Gbs of RAM while guaranteeing a very high performance. We do this by exploiting state-of-the-art bioinformatics tools. DNAscan can screen your DNA NGS data for single nucleotide variants, small indels, structural variants, repeat expansions, viral (or any other organism’s) genetic material. Its results are annotated using a wide range of databases including ClinVar, EXAC, dbSNP and CADD and uploaded onto the gene.iobio platform for an on-the-fly analysis/interpretation

**Version:** 0.1-alpha

## Documentation

[text and links out to ./docs/]  

- description: Aims etc etc
- get the code: git clone   
- dependencies (software, hardware) :<br />
  1)  ExpansionHunter 2.0.9 <br />
  2)  FreeBayes 1.0.2-33<br />
  3)  Manta 1.0.3<br />
  4)  Hisat2 2.1.0<br />
  5)  Python >= 2.7<br />
  6)  Annovar "Version: $Date: 2016-02-01 00:11:18 -0800 (Mon,  1 Feb 2016)"<br />
  7)  Samtools >= 1.3<br />
  8)  Bedtools2 2.25<br />
- test <br />
  python run_analysis.py -format input_file_fomat -reference hg19 -in input_file -out out_dir -SV -BED<br />

## Core Contributors
- [Dr Alfredo Iacoangeli](alfredo.iacoangeli@kcl.ac.uk), UK
- [Dr Stephen J Newhouse](stephen.j.newhouse@gmail.com), UK

### Contributors
- [NAME](email), [COUNTRY OR AFFILIATION]

For a full list of contributors see [LINK](./CONTRIBUTORS.md)

## Contributing

Here’s how we suggest you go about proposing a change to this project:

1. [Fork this project][fork] to your account.
2. [Create a branch][branch] for the change you intend to make.
3. Make your changes to your fork.
4. [Send a pull request][pr] from your fork’s branch to our `master` branch.

Using the web-based interface to make changes is fine too, and will help you
by automatically forking the project and prompting to send a pull request too.

[fork]: https://help.github.com/articles/fork-a-repo/
[branch]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository
[pr]: https://help.github.com/articles/using-pull-requests/


## Licence 
- [MIT](./LICENSE.txt)

## Citation

- [DOI FROM ZENADO]
- [PAPER]

*********



<p align="center">
  Core Developers funded as part of:</br> 
  <a href="https://www.mndassociation.org/">MNDA</a></br> 
  <a href="http://www.maudsleybrc.nihr.ac.uk/">NIHR Maudsley Biomedical Research Centre (BRC), King's College London</a></br>
  <a href="http://www.ucl.ac.uk/health-informatics/">Farr Institute of Health Informatics Research, UCL Institute of Health Informatics, University College London</a>
</p>

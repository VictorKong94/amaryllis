# amaryllis

### Requirements

You should have enough free storage space on your hard drive to store at least
five times the size of your RNA-Seq dataset.
In order to run this project, you must also have the following installed:
- `FastQC`
- `Trimmomatic`
- `Bowtie2`
- `Perl` (modules: `Log::Reproducible`, `autodie.pm`, `Parallel/ForkManager.pm`)
- `R` (modules: `limma`, `edgeR`, `locfit`, `ggplot2`)

If using a server hosted by Amazon Web Services, see
[here](https://gist.github.com/mfcovington/27746b491743ababf32cbadd49846730)
for more information on dependencies.

### Setup and Execution

When configuring the pipeline for the first time, edit the "SETUP" block in
`pipeline.sh` to match your system's directory sturcture and your preferences.

It is highly recommended to perform a test run before using the pipeline on a
full-sized dataset. To do this, run the command `make test-run` in the terminal.

Before running an instance of the pipeline, create a `parameters_*.txt` file
that includes specifications about the job you want to do. See
`parameters_template.txt` for more information.

To run the pipeline, execute the following command in your terminal:
    `sh /path/to/pipeline.sh /path/to/parameters_*.txt`

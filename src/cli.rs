use clap::Parser;


#[derive(Parser, Debug)]
#[command( name = "Diplinator", about = "Diplinator: Choose the best alignment to each haploid of a diploid assembly", version = "0.1.0")]

pub struct Cli {
    // 
    #[arg(short,long, value_name = "FILE", required = true, help="hap1.sam/bam/cram")]
    pub mat: String,

    #[arg(short, long, value_name = "FILE", required = true, help="hap2.sam/bam/cram")]
    pub pat: String, 

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (mat)")]
    pub ref_mat: Option<String>,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (pat)")]
    pub ref_pat: Option<String>,

    #[arg(short, long, value_name = "PREFIX", default_value = "diplinator_out", help="prefix of output files")]
    pub out: String, 

    // inputs are PAF files
    #[arg(long,value_name = "BOOL", default_value_t = false, help = "input files are PAF")]
    pub paf: bool,

    // number of total threads to use
    #[arg(short, long,value_name = "INT", default_value_t = 8, help = "Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance.")]
    pub threads: usize
}
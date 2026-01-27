use clap::Parser;


#[derive(Parser, Debug)]
#[command( name = "diploidinator", about = "Choose the best alignment to each haploid of a diploid assembly", version = "1.1")]

pub struct Cli {
    // 
    #[arg(short,long, value_name = "FILE", required = true, help="hap1.sam/bam/cram")]
    pub mat: String,

    #[arg(short, long, value_name = "FILE", required = true, help="hap2.sam/bam/cram")]
    pub pat: String, 

    // inputs are PAF files
    #[arg(long,value_name = "BOOL", default_value_t = false, help = "input files are PAF")]
    pub paf: bool,

    // number of threads to use
    #[arg(short, long,value_name = "INT", default_value_t = 2, help = "Number of threads to use")]
    pub threads: usize
}
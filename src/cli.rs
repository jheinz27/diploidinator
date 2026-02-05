use clap::Parser;


#[derive(Parser, Debug)]
#[command( name = "diploidinator", about = "Choose the best alignment to each haploid of a diploid assembly", version = "1.1")]

pub struct Cli {
    // 
    #[arg(short,long, value_name = "FILE", required = true, help="hap1.sam/bam/cram")]
    pub mat: String,

    #[arg(short, long, value_name = "FILE", required = true, help="hap2.sam/bam/cram")]
    pub pat: String, 

    //needed for cram format TODO:Add in error handling to make sure this is provifed 
    //NEED REF FOR MAT AND FOR PAT 
    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file")]
    pub ref_mat: Option<String>,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file")]
    pub ref_pat: Option<String>,

    #[arg(short, long, value_name = "PREFIX", default_value = "diploidinator_out", help="hap2.sam/bam/cram")]
    pub out: String, 

    // inputs are PAF files
    #[arg(long,value_name = "BOOL", default_value_t = false, help = "input files are PAF")]
    pub paf: bool,

    // number of threads to use
    // need to divide by 4 and take floor
    #[arg(short, long,value_name = "INT", default_value_t = 4, help = "Number of threads to use for BAM file decompression")]
    pub threads: usize
}
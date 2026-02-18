use clap::Parser;
use diplinator::{Cli, paf, sam};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    let start = Instant::now();
    let args = Cli::parse();

    
    if args.paf {  
         paf::process_paf(&args)?; 

    } else {
        sam::process_sam(&args)?; 
    }

    let duration = start.elapsed();
    eprintln!("Time elapsed: {:?}", duration);
    Ok(())
}


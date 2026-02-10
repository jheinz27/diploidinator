use rust_htslib::{
    bam::{self, Read, Record, Writer, record::Aux },
    errors::Error as BamError,};
use rust_htslib::{bam::record::{Cigar}};
use crate::cli::Cli;
use std::iter::Peekable;
use rust_htslib::{htslib};
use std::ffi::CString;
use std::path::Path;
use std::cmp::max; 


/// Helper function to peek at the file format using c path 
fn get_format_from_path<P: AsRef<Path>>(path: P) -> Result<bam::Format, Box<dyn std::error::Error>> {
    let path_str = path.as_ref().to_str().ok_or("Invalid UTF-8 path")?;
    let c_path = CString::new(path_str).unwrap();

    unsafe {
        let hts_file = htslib::hts_open(c_path.as_ptr(), b"r\0".as_ptr() as *const i8);
        if hts_file.is_null() {
            return Err(format!("Could not open file: {}", path_str).into());
        }
        let format_struct = (*hts_file).format;
        htslib::hts_close(hts_file);

        // Map the C format to the Rust enum
        match format_struct.format {
            htslib::htsExactFormat_bam => Ok(bam::Format::Bam),
            htslib::htsExactFormat_cram => Ok(bam::Format::Cram),
            htslib::htsExactFormat_sam => Ok(bam::Format::Sam), 
            _ => Err(format!("Unsupported or unknown file format for: {}", path_str).into()),
        }
    }

}


fn formats_equal(a: &bam::Format, b: &bam::Format) -> bool {
    match (a, b) {
        (bam::Format::Bam, bam::Format::Bam) => true,
        (bam::Format::Cram, bam::Format::Cram) => true,
        (bam::Format::Sam, bam::Format::Sam) => true,
        _ => false,
    }
}

pub fn process_sam(args: &Cli) -> Result<(), Box<dyn std::error::Error>> { 

    // read in files 
    let mat_format = get_format_from_path(&args.mat).expect("Failed to identify maternal file format");
    let pat_format = get_format_from_path(&args.pat).expect("Failed to identify paternal file format");
    
    if !formats_equal(&mat_format, &pat_format) {
        eprintln!("Error: Input files must have the same format (found {:?} and {:?}).", mat_format, pat_format);
        std::process::exit(1);
    }

    //assert_eq!(mat_format, pat_format, "Input files must have the same format (e.g. both BAM, both SAM, or both CRAM)");

    let mut mat_reader = bam::Reader::from_path(&args.mat).expect("Failed to open -m file");
    let mut pat_reader = bam::Reader::from_path(&args.pat).expect("Failed to open -p file");
    
    let header_mat=  bam::Header::from_template(mat_reader.header());
    let header_pat= bam::Header::from_template(pat_reader.header());

    let extension = match  mat_format{
        bam::Format::Bam => ".bam", 
        bam::Format::Sam => ".sam", 
        bam::Format::Cram => ".cram", 
    };



    let mut out_mat = Writer::from_path (args.out.clone() + "_mat" + extension, &header_mat, mat_format)?;
    let mut out_pat= Writer::from_path (args.out.clone() + "_pat" + extension, &header_pat, pat_format)?;

    if let bam::Format::Cram = mat_format {
        if let Some(reference) = &args.ref_mat { 
            mat_reader.set_reference(reference).expect("Failed to set reference for Maternal Reader"); 
            out_mat.set_reference(reference).expect("Failed to ses reference for Maternal Writer");
        } else {
            eprintln!("Error: input format is CRAM, but no reference FASTA for Maternal Hap provided.");
            eprintln!("Usage hint: add --ref-mat <FILE>");
            std::process::exit(1);
        }
    }

    
    if let bam::Format::Cram = pat_format {
        if let Some(reference) = &args.ref_pat { 
            pat_reader.set_reference(reference).expect("Failed to set reference for Paternal Reader");
            out_pat.set_reference(reference).expect("Failed to set reference for Paternal Writer");
        } else {
            eprintln!("Error: input format is CRAM, but no reference FASTA for Paternal Hap was provided.");
            eprintln!("Usage hint: add --ref-pat <FILE>");
            std::process::exit(1);
        }
    }


    //set threads
    //if user specifies less than 4, set to 4 (1 thread for each reader and each writer is needed)
    let avail_threads = max(4, args.threads); 
    //assign write:reader threads 3:1
    let r = max(1, avail_threads / 8); 
    //if any additional threads available, assign to writer
    //if num threads is odd, leave one idle
    let w = (avail_threads - (2 * r)) / 2;

    //assign threads to each reader/writer pair
    mat_reader.set_threads(r )?;
    pat_reader.set_threads(r)?;
    out_mat.set_threads(w)?;
    out_pat.set_threads(w)?;

    //println!("reader_threads: {}, writer_threads: {}", r, w ); 


    //peakable iterators of each file
    let mut mat_iter = mat_reader.records().peekable();
    let mut pat_iter = pat_reader.records().peekable(); 

    //vectors that store all alignments of one read (cluster of alignments)
    let mut cluster_mat: Vec<Record> = Vec::with_capacity(10); 
    let mut cluster_pat: Vec<Record> = Vec::with_capacity(10);  


    while mat_iter.peek().is_some() || pat_iter.peek().is_some() {
        //move forward by one read for both files
        let _ = get_clusters(&mut mat_iter, &mut cluster_mat);
        let _ = get_clusters(&mut pat_iter, &mut cluster_pat); 
        
        // handle end-of-file / empty cluster cases CHECK
        match (cluster_mat.first(), cluster_pat.first()) {
            (None, None) => break,           // both done
            (Some(_), None) | (None, Some(_)) => {
                return Err("alignment streams out of sync: one file ended earlier".into());
            }
            (Some(m), Some(p)) => {
                if m.qname() != p.qname() {
                    return Err(format!(
                        "alignment streams out of sync: mat={} pat={}",
                        String::from_utf8_lossy(m.qname()),
                        String::from_utf8_lossy(p.qname()),
                    ).into());
                }
            }
        }
        //get cluster wiht the higher alignment score
        let is_winner_mat = compare_clusters(&mut cluster_mat, &mut cluster_pat)?; 
        
        //higher alignment score to maternal hap
        if is_winner_mat {
            for rec in cluster_mat.iter_mut() {
                //output to mat file
                out_mat.write(rec)?;
            }
        } else {
            for rec in cluster_pat.iter_mut() {
                //output to pat file 

                out_pat.write(rec)?;

            }
        }

    }
Ok(())
}

//function to move ahead one read group at a time for SAM/BAM/CRAM
fn get_clusters<I>(records: &mut Peekable<I>, cluster: &mut Vec<Record>)-> Result<(), Box<dyn std::error::Error>>
where 
    I: Iterator<Item= Result<Record,BamError>>,
{
    cluster.clear();

    let first_record = match records.next() {
        Some(Ok(r)) => r,
        Some(Err(e)) => return Err(Box::new(e)), // CRASH HERE if file is bad
        None => return Ok(()), // End of file (Normal)
    };

    let cur_id = first_record.qname().to_vec(); 
    cluster.push(first_record); 

    loop {
        // We peek at next line
        let peek_result = records.peek();

        match peek_result {
            //Next record is valid
            Some(Ok(next_rec)) => {
                if next_rec.qname() == cur_id.as_slice() {
                    // next record belongs to this cluster, consume
    
                    let rec = records.next().unwrap().unwrap();
                    cluster.push(rec);
                } else {
                    // Belongs to the NEXT cluster. 
                    break; 
                }
            },
            // Next record is an ERROR
            Some(Err(_)) => {
            
                let err = records.next().unwrap().unwrap_err();
                return Err(Box::new(err)); 
            },
            //End of File
            None => break,
        }
    }

    return Ok(())
}

//helper function to get weighted score of reads
fn get_weighted_as(cur_clust : &mut Vec<Record>) -> Result<f32, Box<dyn std::error::Error>> {
    
    let mut sum_alignment_lens = 0; 
    let mut sum_alignment_scores = 0;
    let mut read_intervals: Vec<(u32, u32)> = Vec::with_capacity(cur_clust.len());
    
    for rec in cur_clust {
        //do not factor secondary alignments into choosing best alignment
        if rec.is_secondary() {continue}; 

        sum_alignment_lens += get_alignment_len(rec);

        let alignment_score: i32 = match rec.aux(b"AS") {
            Ok(Aux::I8(v))  => v as i32,
            Ok(Aux::I16(v)) => v as i32,
            Ok(Aux::I32(v)) => v,
            Ok(Aux::U8(v))  => v as i32,
            Ok(Aux::U16(v)) => v as i32, 
            Ok(Aux::U32(v)) => v as i32, 
            _ => return Err(format!("Read '{}' is missing the 'AS' tag", 
            String::from_utf8_lossy(rec.qname())).into()),
        };

        
        sum_alignment_scores += alignment_score; 
        let read_start = get_query_start(rec); 
    
        read_intervals.push((read_start, read_start + get_alignment_len(rec))) 
    }
    if sum_alignment_lens == 0 {
        return Ok(0.0);
    } 

    //get total read bases aligned in any record 
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals); 

    return Ok((sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32); 



}


//choose which alignment block to keep 
fn compare_clusters<'a>(clust1:&'a mut Vec<Record>, clust2:&'a mut Vec<Record>) ->  Result<bool,Box<dyn std::error::Error>> {
    
    if clust1.is_empty() || clust2.is_empty() {
        return Err("Fatal Error: Attempted to compare empty read clusters. This usually indicates a file sync issue.".into());
    }

    let unmappeds = (clust1[0].is_unmapped(), clust2[0].is_unmapped()); 

    match unmappeds {
        (true, true) => {
            let qname_string = String::from_utf8_lossy(clust1[0].qname()).into_owned();
            return Ok(crate::choose_random(qname_string));
        }, //both unmapped 
        (true, false) => return Ok(false), // pat mapped 
        (false, true) => return Ok(true), // mat mapped
        _ => {} //both mapped
    }

    //get AS score for read to each haploid

    let score1 = get_weighted_as(clust1)?; 
    let score2 = get_weighted_as(clust2)?; 

    //return bool of higher alignment
    if score1 > score2 {
        Ok(true)
    } else if score1 < score2 {
        Ok(false)
    } else {
        let qname_string = String::from_utf8_lossy(clust1[0].qname()).into_owned();
        Ok(crate::choose_random(qname_string))
    }
} 


//function to get query span of aligned seqment
fn get_alignment_len(rec: &Record) -> u32  {
    let mut qlen = 0;
    //parse cigar string to determine total ALIGNED query length 
    for c in rec.cigar().iter() {
        match c {
            //consumes query  
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l) => {qlen += *l},
            _ => {}
        }
    }
    return qlen; 
}

//function to get start of query span 
fn get_query_start(rec: &Record) -> u32 {
    let cigar = rec.cigar();
    if rec.is_reverse() {
        let mut right = 0;
        for c in cigar.iter().rev() {
            match *c {
                Cigar::HardClip(l) | Cigar::SoftClip(l) => right += l,
                _ => break, 
            }
        }
        return right; 
    } else {
        let mut left = 0;
        for c in cigar.iter() {
            match *c {
                Cigar::HardClip(l) | Cigar::SoftClip(l) => left += l,
                _ => break,
            }
        }
        return left; 
    }
    

}


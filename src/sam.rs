use rust_htslib::{
    bam::{self, Read, Record, Writer, Format, record::Aux },
    errors::Error as BamError,};
use rust_htslib::{bam::record::{Cigar}};
use crate::cli::Cli;
use std::iter::Peekable;



pub fn process_sam(args: &Cli) -> Result<(), Box<dyn std::error::Error>> { 

    // read in files 
    let mut mat_reader = bam::Reader::from_path(&args.mat).expect("Failed to open -m file");
    let mut pat_reader = bam::Reader::from_path(&args.pat).expect("Failed to open -p file");


    let header_mat=  bam::Header::from_template(mat_reader.header());
    let header_pat= bam::Header::from_template(pat_reader.header());
    let mut out_mat = Writer::from_path ("mat_reads.sam", &header_mat, Format::Sam)?;
    let mut out_pat= Writer::from_path ("pat_reads.sam", &header_pat, Format::Sam)?;

    //set threads
    mat_reader.set_threads(args.threads)?;
    pat_reader.set_threads(args.threads)?;
    //out_shared.set_threads(args.threads)?;

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
                return Err("SAM streams out of sync: one file ended earlier".into());
            }
            (Some(m), Some(p)) => {
                if m.qname() != p.qname() {
                    return Err(format!(
                        "SAM streams out of sync: mat={} pat={}",
                        String::from_utf8_lossy(m.qname()),
                        String::from_utf8_lossy(p.qname()),
                    ).into());
                }
            }
        }
        //get cluster wiht the higher alignment score
        let is_winner_mat = compare_clusters(&mut cluster_mat, &mut cluster_pat); 
        
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
fn get_clusters<I>(records: &mut Peekable<I>, cluster: &mut Vec<Record>)-> Result<(), BamError>
where 
    I: Iterator<Item= Result<Record,BamError>>,
{
    cluster.clear();
    if let Some(Ok(record)) = records.next() {
        let cur_id = record.qname().to_vec(); 
        cluster.push(record); 
        
        while let Some(Ok(next)) = records.peek() {
            if  next.qname() == cur_id.as_slice() {
                if let Some(rec)= records.next() {
                    cluster.push(rec?);
                }
            } else {
                break;
            }
        }
    };
    return Ok(())
}

//helper function to get weighted score of reads
fn get_weighted_as(cur_clust : &mut Vec<Record>) -> f32 {
    
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
            _ => 0,// TODO: FIX to throw an error here
        };

        
        sum_alignment_scores += alignment_score; 
        let read_start = get_query_start(rec); 
        //TODO: check handling of reverse strand cases, may need to change alignment_len
        read_intervals.push((read_start, read_start + get_alignment_len(rec))) 
    }
    //get total read bases aligned in any record 
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals); 
    println!("al: {}", sum_alignment_lens); 
    println!("as: {}", sum_alignment_scores); 
    println!("intervals: {:?}", read_intervals); 

    //weighted_score = (SUM(Alignment_Score) / SUM(Alignment_len)) * tot read_bps_aligned 
    return (sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32; 



}


//choose which alignment block to keep 
fn compare_clusters<'a>(clust1:&'a mut Vec<Record>, clust2:&'a mut Vec<Record>) ->  bool {
    
    //TODO: add error handling for empty clusters
    let unmappeds = (clust1[0].is_unmapped(), clust2[0].is_unmapped()); 

    match unmappeds {
        (true, true) => {
            let qname_string = String::from_utf8_lossy(clust1[0].qname()).into_owned();
            return crate::choose_random(qname_string);
        }, //both unmapped 
        (true, false) => return false, // pat mapped 
        (false, true) => return true, // mat mapped
        _ => {} //both mapped
    }

    //get AS score for read to each haploid
    println!("hap1"); 
    let score1 = get_weighted_as(clust1); 
    println!("s1:{}", score1); 
    println!("hap2"); 
    let score2 = get_weighted_as(clust2); 
    println!("s2:{}", score2); 

    //return bool of higher alignment
    if score1 > score2 {
        true
    } else if score1 < score2 {
        false
    } else {
        let qname_string = String::from_utf8_lossy(clust1[0].qname()).into_owned();
        crate::choose_random(qname_string)
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


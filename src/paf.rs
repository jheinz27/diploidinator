use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write, BufWriter},
    iter::Peekable,
};
use crate::cli::Cli;

//TODO: handle secondary alingments!!! 

pub fn process_paf(args: &Cli) -> Result<(), Box<dyn std::error::Error>> {

    //let stdout = io::stdout();
    //let mut handle = stdout.lock();

    let file1 = File::open( &args.mat)?;
    let file2 = File::open(&args.pat)?;

    let mut reader1 = BufReader::new(file1).lines().peekable();
    let mut reader2 = BufReader::new(file2).lines().peekable();

    let mut out_mat = BufWriter::new(File::create("out_mat.paf")?);
    let mut out_pat = BufWriter::new(File::create("out_pat.paf")?);

    let mut cluster_mat = Vec::with_capacity(10); 
    let mut cluster_pat = Vec::with_capacity(10); 

    while let Some(Ok(_)) = reader1.peek() {
        let _ = get_clusters(&mut reader1, &mut cluster_mat);
        let _ = get_clusters(&mut reader2, &mut cluster_pat);

        //get cluster wiht the higher weighted alignment score
        let is_winner_mat = compare_clusters(&mut cluster_mat, &mut cluster_pat); 
        
        // write all records of read with higher weighted alignment score
        if is_winner_mat {
            for rec in cluster_mat.iter_mut() {
                //output to mat file
                writeln!(out_mat, "{}", rec)?;
            }
        } else {
            for rec in cluster_pat.iter_mut() {
                //output to pat file 
                writeln!(out_pat, "{}", rec)?;

            }
        }

    }
    Ok(())
}

//function to move ahead one read group at a time
pub fn get_clusters<I>(lines: &mut Peekable<I>, cluster: &mut Vec<String>)-> io::Result<()>
where 
    I: Iterator<Item= Result<String, std::io::Error>>,
{
    cluster.clear();
    if let Some(Ok(line)) = lines.next() {
        let cur_id = line.split_once('\t').unwrap().0.to_string(); 
        cluster.push(line); 
        
        while let Some(Ok(next)) = lines.peek() {
            let next_id = next.split_once('\t').unwrap().0;
            if  next_id == cur_id {
                cluster.push(lines.next().unwrap()?)
            } else {
                break;
            }
        }
    };
    return Ok(())
}



//helper function to get weighted score of split reads
//weighted_score = (SUM(AS) / SUM(Alingment len)) * tot read_bps_aligned 
pub fn get_weighted_as(cur_clust : &Vec<String>) -> f32 {
    let mut sum_alignment_lens = 0; 
    let mut sum_alignment_scores = 0;
    let mut read_intervals: Vec<(u32, u32)> = Vec::with_capacity(cur_clust.len());

    for alignment in cur_clust {
        let fields: Vec<&str> = alignment.split('\t').collect(); 
        
        //do not factor secondary alignments into choosing best alignment
        let is_secondary = fields.iter()
        .find(|&&f| f.starts_with("tp:A:"))
        .map(|s| s.ends_with('S')) 
        .unwrap_or(false);

        if is_secondary { continue; } 

        //get alignment length in read_coordinates
        let read_start: u32 = fields[2].parse().unwrap(); 
        let read_end: u32 = fields[3].parse().unwrap(); 
        sum_alignment_lens += read_end - read_start;

        //store read alignment coordinates to merge all overlaps at end
        read_intervals.push((read_start, read_end));

        //extract alignment score in AS tag and add to running sum 
        let as_score = fields.iter()
            .find(|&&f| f.starts_with("AS:i:"))
            .map(|s| s[5..].parse::<i32>().unwrap_or(0))
            .unwrap_or(0); // default to 0 if tag missing 
        sum_alignment_scores += as_score; 

    }
     
    //get total read bases aligned in any record 
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals); 

    //weighted_score = (SUM(Alignment_Score) / SUM(Alignment_len)) * tot read_bps_aligned 
    return (sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32; 

}


pub fn compare_clusters<'a>(clust1:&'a Vec<String>, clust2:&'a Vec<String>) -> bool{
    match (clust1[0].split('\t').nth(4), clust2[0].split('\t').nth(4)) {
        (Some("*"), Some("*")) => return crate::choose_random(clust1[0].split('\t').next().unwrap().to_string()), // both reads unmapped
        (Some("*"), _) => return false, // mat hap unmapped
        (_, Some("*")) => return true, // pat ap unmapped
        _ => {} // Continue if mapped to both haps 
    }

    let score1 = get_weighted_as(&clust1); 
    let score2 = get_weighted_as(&clust2); 

    //return bool of higher alignment
    if score1 > score2 {
        true
    } else if score1 < score2 {
        false
    } else {
        crate::choose_random(clust1[0].split('\t').next().unwrap().to_string())
    }
} 
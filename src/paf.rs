use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write, BufWriter},
    iter::Peekable,
};
use crate::cli::Cli;


pub fn process_paf(args: &Cli) -> Result<(), Box<dyn std::error::Error>> {

    //let stdout = io::stdout();
    //let mut handle = stdout.lock();

    let file1 = File::open( &args.mat)?;
    let file2 = File::open(&args.pat)?;

    let mut mat_iter = BufReader::new(file1).lines().peekable();
    let mut pat_iter = BufReader::new(file2).lines().peekable();

    
    let mut out_mat = BufWriter::new(File::create( args.out.clone() + "_mat.paf")?);
    let mut out_pat = BufWriter::new(File::create( args.out.clone()+ "_pat.paf")?);

    let mut cluster_mat = Vec::with_capacity(10); 
    let mut cluster_pat = Vec::with_capacity(10); 

    while mat_iter.peek().is_some() || pat_iter.peek().is_some() {
        let _ = get_clusters(&mut mat_iter, &mut cluster_mat);
        let _ = get_clusters(&mut pat_iter, &mut cluster_pat);

        let id1 = cluster_mat[0].split('\t').next().unwrap();
        let id2 = cluster_pat[0].split('\t').next().unwrap();
        assert_eq!(id1, id2, "PAF streams out of sync: {id1} vs {id2}");
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
        let mut fields = alignment.split('\t');

        //skip to relevant columns
        let _qname = fields.next();
        let _qlen = fields.next();
        let qstart = fields.next().unwrap().parse::<u32>().unwrap();
        let qend = fields.next().unwrap().parse::<u32>().unwrap();

        let mut is_secondary = false; 
        let mut as_score = 0;
        //find AS and tp tags (start from field 12)
        for field in fields.skip(8) { // skip to tags
            if field.starts_with("tp:A:S") { is_secondary = true; }
            if field.starts_with("AS:i:") { 
                 as_score = field[5..].parse().unwrap_or(0);
            }
        }
        
        //do not factor secondary alignments into choosing best alignment
        if is_secondary { continue; } 
    

        sum_alignment_lens += qend - qstart;

        //store read alignment coordinates to merge all overlaps at end
        read_intervals.push((qstart, qend));

        sum_alignment_scores += as_score; 

    }
     
    //get total read bases aligned in any record 
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals); 

    //weighted_score = (SUM(Alignment_Score) / SUM(Alignment_len)) * tot read_bps_aligned 
    return (sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32; 

}


pub fn compare_clusters<'a>(clust1:&'a Vec<String>, clust2:&'a Vec<String>) -> bool{
    match (clust1[0].split('\t').nth(5), clust2[0].split('\t').nth(5)) {
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
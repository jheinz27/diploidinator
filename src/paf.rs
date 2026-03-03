use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write, BufWriter},
    iter::Peekable,
};
use crate::cli::Cli;


pub fn process_paf(args: &Cli) -> Result<(), Box<dyn std::error::Error>> {

    // read in both files
    let file1 = File::open(&args.asm1)
        .map_err(|e| format!("Failed to open asm1 file '{}': {}", args.asm1, e))?;
    let file2 = File::open(&args.asm2)
        .map_err(|e| format!("Failed to open asm2 file '{}': {}", args.asm2, e))?;

    //create peekable iterators of each file (line-by-line for PAF)
    let mut asm1_iter = BufReader::new(file1).lines().peekable();
    let mut asm2_iter = BufReader::new(file2).lines().peekable();

    //create writers for both outputs that share user specified prefix
    let asm1_out_path = format!("diplinator_{}.paf", args.s1);
    let asm2_out_path = format!("diplinator_{}.paf", args.s2);
    let mut out_asm1 = BufWriter::new(File::create(&asm1_out_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", asm1_out_path, e))?);
    let mut out_asm2 = BufWriter::new(File::create(&asm2_out_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", asm2_out_path, e))?);

    //vectors that store all alignments of one read (cluster of alignments)
    //initialize capacity to 10 to account for supplemental and secondary alignments
    let mut cluster_asm1: Vec<String> = Vec::with_capacity(10);
    let mut cluster_asm2: Vec<String> = Vec::with_capacity(10);

    //initialize counts for summary statistics printed to terminal
    let mut count_asm1: u64 = 0;
    let mut count_asm2: u64 = 0;
    let mut count_equal: u64 = 0;
    let mut count_unmapped: u64 = 0;

    //iterate through both files until they are both exhausted
    while asm1_iter.peek().is_some() || asm2_iter.peek().is_some() {

        //move forward by one read for both files
        get_clusters(&mut asm1_iter, &mut cluster_asm1)?;
        get_clusters(&mut asm2_iter, &mut cluster_asm2)?;

        // check for possible errors such as:
        //end of file / empty cluster / clusters don't represent same read in both files
        match (cluster_asm1.first(), cluster_asm2.first()) {
            (None, None) => break,           // end of file reached for both, should occur at same iteration
            (Some(_), None) | (None, Some(_)) => {
                //one file has ended earlier than the other- throw error
                return Err("PAF streams out of sync: one file ended earlier".into());
            }
            (Some(m), Some(p)) => {
                //read ID is not the same in both clusters- throw error
                let id1 = m.split('\t').next().unwrap_or("");
                let id2 = p.split('\t').next().unwrap_or("");
                if id1 != id2 {
                    return Err(format!(
                        "PAF streams out of sync: asm1={} asm2={}", id1, id2
                    ).into());
                }
            }
        }

        //get cluster with the higher alignment score, returns the Winner enum
        let winner = compare_clusters(&mut cluster_asm1, &mut cluster_asm2, args)?;

        //logic for which file to write read to given score comparison output
        match winner {
            //asm1 clear winner, write to out_asm1
            crate::Winner::Asm1 => {
                count_asm1 += 1; // increment read counter
                for rec in cluster_asm1.iter_mut() {
                    writeln!(out_asm1, "{}", rec)?;
                }
            }
            //asm2 clear winner, write to out_asm2
            crate::Winner::Asm2 => {
                count_asm2 += 1; // increment read counter
                for rec in cluster_asm2.iter_mut() {
                    writeln!(out_asm2, "{}", rec)?;
                }
            }
            crate::Winner::Both => {
                count_equal += 1; // increment read counter
                //if user specifies --both, write equal scoring reads to both output files
                if args.both {
                    for rec in cluster_asm1.iter_mut() {
                        writeln!(out_asm1, "{}", rec)?;
                    }
                    for rec in cluster_asm2.iter_mut() {
                        writeln!(out_asm2, "{}", rec)?;
                    }
                //default behavior is randomly assign equal scoring read to one file
                } else {
                    //hash read name and use last bit value to assign to asm1 or asm2
                    //ensures that assignments will be reproducible
                    let qname = cluster_asm1[0].split('\t').next().unwrap();
                    match crate::choose_random(qname.as_bytes()) {
                        crate::Winner::Asm1 => {
                            for rec in cluster_asm1.iter_mut() {
                                writeln!(out_asm1, "{}", rec)?;
                            }
                        }
                        _ => {
                            for rec in cluster_asm2.iter_mut() {
                                writeln!(out_asm2, "{}", rec)?;
                            }
                        }
                    }
                }
            }
            crate::Winner::Unmapped => {
                count_unmapped += 1;
                match args.unmapped {
                    crate::cli::UnmappedDest::Asm1 => {
                        for rec in cluster_asm1.iter_mut() { writeln!(out_asm1, "{}", rec)?; }
                    }
                    crate::cli::UnmappedDest::Asm2 => {
                        for rec in cluster_asm2.iter_mut() { writeln!(out_asm2, "{}", rec)?; }
                    }
                    crate::cli::UnmappedDest::Discard => {}
                }
            }
        }

    }
    //print summary statistics to terminal
    let total = count_asm1 + count_asm2 + count_equal + count_unmapped;
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s1, count_asm1, count_asm1 as f64 / total as f64 * 100.0);
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s2, count_asm2, count_asm2 as f64 / total as f64 * 100.0);
    eprintln!("Reads with equal scores:     {} ({:.1}%)", count_equal, count_equal as f64 / total as f64 * 100.0);
    eprintln!("Reads unmapped to both:      {} ({:.1}%)", count_unmapped, count_unmapped as f64 / total as f64 * 100.0);
    eprintln!("Total reads parsed:          {}", total);
    Ok(())
}

//function to move ahead one read group at a time for PAF
fn get_clusters<I>(lines: &mut Peekable<I>, cluster: &mut Vec<String>)-> io::Result<()>
where
    I: Iterator<Item= Result<String, std::io::Error>>,
{
    //forget previous cluster
    cluster.clear();

    //access alignment record of next line in iterator if it exists
    let first_line = match lines.next() {
        Some(Ok(line)) => line,
        Some(Err(e)) => return Err(e), //throw error if file appears corrupted
        None => return Ok(()),         // End of file
    };

    //get read ID of record (first tab-delimited field in PAF)
    //we cluster any records with the same read ID
    let cur_id = first_line.split_once('\t')
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Malformed PAF line: no tab delimiter"))?.0.to_string();
    //store first record
    cluster.push(first_line);

    //look for further lines with same read ID
    loop {
        //peek at next line
        match lines.peek() {
            //Next record is valid
            Some(Ok(next)) => {
                //check if next record has same read ID
                let next_id = match next.split_once('\t') {
                    Some((id, _)) => id,
                    None => break, // malformed line, let caller handle it
                };
                if next_id == cur_id {
                    // next record belongs to this cluster, consume and add to cluster
                    let line = lines.next().unwrap()?;
                    cluster.push(line);
                } else {
                    // Belongs to the next cluster.
                    break;
                }
            },
            // Next record is corrupt
            Some(Err(_)) => {
                let err = lines.next().unwrap().unwrap_err();
                return Err(err);
            },
            //end of file
            None => break,
        }
    }

    //mutated cluster vector in place, only need to return result Ok
    return Ok(())
}


//helper function to get weighted score of split reads using a specified tag (AS or ms)
//weighted_score = (SUM(score) / SUM(Alignment_len)) * tot read_bps_aligned
pub fn get_weighted_score(cur_clust : &Vec<String>, tag_prefix: &str) -> Result<f32, Box<dyn std::error::Error>> {
    let mut sum_alignment_lens = 0;
    let mut sum_alignment_scores = 0;
    let mut read_intervals: Vec<(u32, u32)> = Vec::with_capacity(cur_clust.len());

    for alignment in cur_clust {
        let mut fields = alignment.split('\t');

        //skip to relevant columns
        let _qname = fields.next().ok_or("Malformed PAF line: missing qname field")?;
        let _qlen = fields.next().ok_or("Malformed PAF line: missing qlen field")?;
        let qstart = fields.next()
            .ok_or("Malformed PAF line: missing qstart field")?
            .parse::<u32>().map_err(|e| format!("Invalid qstart in PAF: {}", e))?;
        let qend = fields.next()
            .ok_or("Malformed PAF line: missing qend field")?
            .parse::<u32>().map_err(|e| format!("Invalid qend in PAF: {}", e))?;

        let mut is_secondary = false;
        let mut as_score = 0;
        //find score and tp tags (start from field 12)
        for field in fields.skip(8) { // skip to tags
            if field.starts_with("tp:A:S") { is_secondary = true; } // check for secondary alignment tag
            if field.starts_with(tag_prefix) {
                 as_score = field[tag_prefix.len()..].parse()
                    .map_err(|e| format!("Invalid {} value in PAF: {}", tag_prefix, e))?;
            }
        }

        //do not factor secondary alignments into choosing best alignment
        if is_secondary { continue; }

        //should not happen if paf is formatted correctly
        if qend < qstart {
            return Err(format!("Malformed PAF: qend ({}) < qstart ({}) for read '{}'",
                qend, qstart, cur_clust[0].split('\t').next().unwrap_or("unknown")).into());
        }

        sum_alignment_lens += qend - qstart;

        //store read alignment coordinates to merge all overlaps at end
        read_intervals.push((qstart, qend));

        sum_alignment_scores += as_score;

    }

    if sum_alignment_lens <= 0 {
        return Err(format!("Read '{}' has primary alignment length of 0",
            cur_clust[0].split('\t').next().unwrap_or("unknown")).into());
    }

    //get total read bases aligned in any record
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals);

    //weighted_score = (SUM(Alignment_Score) / SUM(Alignment_len)) * tot read_bps_aligned
    return Ok((sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32);

}

pub fn compare_clusters<'a>(clust1:&'a Vec<String>, clust2:&'a Vec<String>, args: &Cli) ->  Result<crate::Winner, Box<dyn std::error::Error>> {

    match (clust1[0].split('\t').nth(5), clust2[0].split('\t').nth(5)) {
        (Some("*"), Some("*")) => {return Ok(crate::Winner::Unmapped);}, // both reads unmapped
        (Some("*"), _) => return Ok(crate::Winner::Asm2), // asm1 hap unmapped
        (_, Some("*")) => return Ok(crate::Winner::Asm1), // asm2 hap unmapped
        _ => {} // continue if mapped to both haps
    }

    let tag_prefix = if args.ms { "ms:i:" } else { "AS:i:" };
    let (score1, score2) = (
        get_weighted_score(clust1, tag_prefix)?,
        get_weighted_score(clust2, tag_prefix)?,
    );

    //return respective winner depending on which AS is higher,
    //both is a special case that can be determined by user input
    if score1 > score2 {
        Ok(crate::Winner::Asm1)
    } else if score1 < score2 {
        Ok(crate::Winner::Asm2)
    } else {
        Ok(crate::Winner::Both)
    }
}

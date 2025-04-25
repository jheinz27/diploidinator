use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::iter::Peekable;
use std::time::Instant;
use std::fs::OpenOptions;

//function to move ahead one read group at a time
fn get_clusters<I>(lines: &mut Peekable<I>, cluster: &mut Vec<String>)-> io::Result<()>
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

//helper function to get score of split reads
//TODO fix to be same as python version (multiplication by length of reads)
fn get_as(cur_clust : &Vec<String>) -> i32 {
    let mut sum = 0; 
    for alignment in cur_clust {
        let field = alignment.split('\t').nth(14).unwrap(); 
        sum += field.split(':').nth(2).unwrap().parse::<i32>().unwrap(); 
    }
    return sum 
}

fn compare_clusters<'a>(clust1:&'a Vec<String>, clust2:&'a Vec<String>) -> Option<&'a Vec<String>> {
    match (clust1[0].split('\t').nth(4), clust2[0].split('\t').nth(4)) {
        (Some("*"), Some("*")) => return None, // both reads unmapped
        (Some("*"), _) => return Some(clust2), // hap1 unmapped
        (_, Some("*")) => return Some(clust1), // hap2 unmapped
        _ => {} // Continue if mapped to both haps 
    }
    
    let score1 = get_as(&clust1); 
    let score2 = get_as(&clust2); 

    Some(if score1 >= score2 {clust1} else {clust2})
} 

fn main() -> io::Result<()> {
    let start = Instant::now(); 
    let args: Vec<String> = env::args().collect();
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let file1 = File::open(&args[1])?;
    let file2 = File::open(&args[2])?;

    let mut reader1 = BufReader::new(file1).lines().peekable();
    let mut reader2 = BufReader::new(file2).lines().peekable();

    let mut cluster1 = Vec::with_capacity(10); 
    let mut cluster2 = Vec::with_capacity(10); 

    while let Some(Ok(_)) = reader1.peek() {
        let _ = get_clusters(&mut reader1, &mut cluster1);
        let _ = get_clusters(&mut reader2, &mut cluster2);
        if let Some(winner_clust) = compare_clusters(&cluster1, &cluster2) {
            writeln!(handle, "{}", winner_clust.join("\n"))?;
        }
    }
    let mut file = OpenOptions::new()
    .create(true)
    .append(true)
    .open("rust_timing.log")?;
    let duration = start.elapsed();
    writeln!(file, "{}\t{:?}", &args[1], duration)?; 
    Ok(())
}
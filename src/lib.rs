pub mod cli;
pub use cli::Cli;
pub mod paf;
pub mod sam;
use std::hash::{Hash, Hasher};
use twox_hash::XxHash64;


//if read has identical alignment to both haps, chose which hap to report randomly
//based on last digit of hash of read ID
//XxHash64 provides reproducible assignment
pub fn choose_random(id: String) -> bool {
    let mut hasher = XxHash64::with_seed(42);
    id.hash(&mut hasher); 
    return hasher.finish() & 1 == 0 ; 
}

//helper function to merge any read alignment segments that overlap in read coordinates
//returns count of unique bps of the read contained in an alignment segment
pub fn merge_intervals(intervals: &mut Vec<(u32, u32)>) -> u32 {
    //sort cluster by read start location of alignment segment
    intervals.sort_unstable_by_key(|k| k.0);
    
    let mut read_bps_aligned = 0; 
    if !intervals.is_empty() {
        //initialize at first interval
        let (mut cur_start, mut cur_end) = intervals[0]; 
        //iterate through intervals and merge adjacent overlapping intervals
        for &(next_start, next_end) in intervals.iter().skip(1) { 
            if next_start < cur_end {
                // intervals overlap, so extend
                if next_end > cur_end {
                    cur_end = next_end;
                }
            } else {
                //no further overlap, add length and start over with next interval grouping
                read_bps_aligned += cur_end - cur_start;
                cur_start = next_start;
                cur_end = next_end;
            }
        
        }
        read_bps_aligned += cur_end - cur_start
    }
    return read_bps_aligned; 
}
 
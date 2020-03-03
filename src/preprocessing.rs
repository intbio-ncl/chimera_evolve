use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use suffix::SuffixTable;

#[allow(dead_code)]
// Given a hash map and codon string, converts the codon string t
pub fn translate_codon_str_to_alphabet(codon_string: &str, map: &HashMap<&str, &str>) -> String {
    codon_string
        .chars()
        .enumerate()
        .filter(|(idx, _)| idx % 3 == 0)
        .map(|(idx, _)| *map.get(&codon_string[idx..idx + 3]).unwrap())
        .collect::<String>()
}

fn condensed_encoding_substring(
    encoded_substring: &SuffixTable,
    aa_seq: &str,
    map: &HashMap<&str, Vec<&str>>,
) -> HashSet<String> {
    // This function gets the relevant substrings encoding the start of a given substring.
    // Create an iterator for the substring.
    let mut aa_iter = aa_seq.chars();

    let start = aa_iter.next();
    // Check to ensure this is not None.
    match start {
        Some(_t) => {}
        None => panic!("The condensed encoding substring should not be empty."),
    }

    // Get an initial vector
    let initial_vec: &Vec<&str> = map.get::<str>(&start.unwrap().to_string()).unwrap();
    // Make a new vector, that contains the strings to expand.
    let mut current_vec: Vec<String> = c![i.to_string(), for i in initial_vec];
    // Make a set to hold all possible encoding substrings.
    let mut all_substrings: HashSet<String> = HashSet::new();

    loop {
        // Update current vector.
        current_vec = current_vec
            .iter()
            // Filter out members that are not in the encoded substring.
            .filter(|i| encoded_substring.contains(*i))
            // Add remaining substrings to all_substrings.
            .map(|i| {
                all_substrings.insert(i.to_string());
                i.to_string()
            })
            // And collect them in a vector
            .collect::<Vec<String>>();


        // If the current vector is broken (nothing to extend), then break.
        if current_vec.is_empty() {
            break;
        }
        // Get the next char
        let next_char = aa_iter.next();
        // If there's no next char, break the loop.
        match next_char {
            Some(_t) => {}
            None => break,
        }

        // Update current_vec to create new possible strings up to the next amino acid.
        current_vec = current_vec
            .iter()
            .map(|i| {
                let mut tmp_vec: Vec<String> = vec![];
                let options = map.get::<str>(&next_char.unwrap().to_string()).unwrap();
                for item in options.iter() {
                    let new_item = format!("{}{}", i, item);
                    tmp_vec.push(new_item);
                }
                tmp_vec
            })
            .flatten()
            .collect::<Vec<String>>();
    }

    all_substrings
        .iter()
        .filter(|i| {
            for j in &all_substrings {
                if j.contains(*i) && j != *i {
                    return false;
                }
            }
            true
        })
        .cloned()
        .collect::<HashSet<String>>()
}

pub fn condense_encoding<'a, 'b>(
    encoded_cds: &str,
    aa_seq: &str,
    map: &HashMap<&str, Vec<&str>>,
) -> SuffixTable<'a, 'b> {
    // This function creates a condensed encoding based on the encoded CDS and AA seq.
    // Create a suffix table
    info!("Generating suffix table for full CDS");
    let st = SuffixTable::new(encoded_cds);

    // Get all the relevant substrings
    info!("Obtaining relevant substrings from host CDS, given protein");
    let mut substrings = aa_seq
        .chars()
        .enumerate()
        .par_bridge()
        .map(|(idx, _)| condensed_encoding_substring(&st, &aa_seq[idx..], &map))
        .flatten()
        .collect::<Vec<String>>();
    info!("Compressing substrings into single string");
    substrings.sort_by(|a, b| a.len().cmp(&b.len()));

    let mut final_string = String::new();
    for string in substrings.iter().rev() {
        if !final_string.contains(string) {
            final_string.push_str(string);
            final_string.push('|');
        }
    }
    info!("Final string generated");
    info!("Creating suffix table");
    SuffixTable::new(final_string)
}

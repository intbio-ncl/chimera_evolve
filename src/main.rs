#[macro_use]
extern crate log;
extern crate env_logger;
extern crate suffix;
#[macro_use(c)]
extern crate cute;
extern crate bio;
extern crate rand;
extern crate rayon;
#[macro_use]
extern crate clap;

use clap::App;
use std::collections::HashMap;
use std::fs::File;
use std::fmt::Error;
use std::io::{BufWriter, Write};
use std::time::{Duration, Instant};


mod ealgorithm;
mod preprocessing;
mod translation_tables;

fn parse_protein(path: &str) -> String {
    let reader = bio::io::fasta::Reader::from_file(path).unwrap();
    let mut sequence = String::new();
    for (idx, item) in reader.records().enumerate() {
        let result = item.unwrap();
        if idx > 0 {
            panic!("Too many records in the given FASTA file for protein.");
        }
        sequence.push_str(std::str::from_utf8(result.seq()).unwrap());
    }
    sequence
}

fn parse_cds(path: &str, codon_singular_map: &HashMap<&str, &str>) -> String {
    info!("{}", format!("Parsing: {}", path));
    let reader = bio::io::fasta::Reader::from_file(path).unwrap();
    let mut encoding = String::new();
    for item in reader.records() {
        let result = item.unwrap();
        let seq = std::str::from_utf8(result.seq()).unwrap();
        let encoded_seq = preprocessing::translate_codon_str_to_alphabet(seq, &codon_singular_map);
        encoding.push_str(&encoded_seq);
        encoding.push('|');
    }

    encoding
}

fn setup_logger() {
    std::env::set_var("RUST_LOG", "INFO");
    env_logger::init();
}

fn encoding_to_nuc(encoded_str: &str, map: &HashMap<&str, &str> ) -> String {
    let mut nucleotide_string = String::new();
    encoded_str.chars().for_each(|i| {
        let i_str = i.to_string();
        nucleotide_string.push_str( map.get::<str>(&i_str).unwrap() );
    });
    nucleotide_string
}

fn write_to_fasta(outfile: &str, sequence: &str, score: f64, duration: &Duration) -> Result<(), Error> {
    // Open the outfile and create a buffer.
    let file = File::create(outfile).unwrap();
    let mut buf = BufWriter::new(file);
    let header = format!(">Result score={} duration={:?}\n", score, duration);

    let mut check = buf.write(header.as_bytes());
    check.unwrap();

    for (pos, item) in sequence.chars().enumerate() {

        if pos % 60 == 0 && pos != 0 {
            check = buf.write("\n".as_bytes());
            check.unwrap();
        }

        check = buf.write(item.to_string().as_bytes());
        check.unwrap();


    }
    buf.flush().unwrap();

    Ok(())

}

fn main() {
    // Parse command line arguments
    let yaml = load_yaml!("cli.yml");
    let matches = App::from_yaml(yaml).get_matches();

    info!("Parsing command line arguments");
    // Main files.
    let genomes: Vec<_> = matches.values_of("cds").unwrap().collect();
    let protein: &str = matches.value_of("protein").unwrap();
    // Algorithm parameters (numeric)
    let no_mutations: usize = matches
        .value_of("mutations")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let no_crossovers: usize = matches
        .value_of("crossovers")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let no_generations: usize = matches
        .value_of("generations")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let gen_start: usize = matches
        .value_of("generation_start")
        .unwrap()
        .parse::<usize>()
        .unwrap();
    let outfile: &str = matches
        .value_of("outfile")
        .unwrap();
    // Method and weights.
    let method = matches.value_of("method").unwrap();
    let _weights = matches.value_of("weights");
    let weights: Option<Vec<f64>> = match _weights {
        Some(t) => Some(
            t.split(',')
                .filter_map(|s| s.parse::<f64>().ok())
                .collect::<Vec<_>>(),
        ),
        None => None,
    };

    // Set up the logger.
    setup_logger();
    info!("Script started");

    let start = Instant::now();

    info!("Generating maps for codon table 11");
    let (
        codon_singular_map,
        singular_codon_map,
        _codon_amino_map,
        _amino_codon_map,
        amino_singular_map,
        alternatives,
    ) = translation_tables::generate_maps(translation_tables::tt11());

    info!("Parsing the protein sequence");
    let pro_seq = parse_protein(protein);
    info!("Parsing and encoding coding sequences");
    let suffix_tables = c![
        preprocessing::condense_encoding(
            &parse_cds(i, &codon_singular_map), &pro_seq, &amino_singular_map
        ), for i in genomes];

    // Run the genetic algorithm&
    let (result_seq, result_fitness) = ealgorithm::run_ea(
        &pro_seq,
        &suffix_tables,
        &amino_singular_map,
        &alternatives,
        no_crossovers,
        no_mutations,
        no_generations,
        gen_start,
        method,
        &weights,
    );

    let encoded = encoding_to_nuc(&result_seq, &singular_codon_map);
    let duration = start.elapsed();
    write_to_fasta(outfile, &encoded, result_fitness, &duration).unwrap();
    info!("Chimera evolve algorithm completed -- thank you for flying Air ICOS! 🚀")
}

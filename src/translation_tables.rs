use std::collections::HashMap;

type StrStrMap = HashMap<&'static str, &'static str>;
type StrStrVecMap = HashMap<&'static str, Vec<&'static str>>;

#[allow(dead_code)]
pub fn tt11() -> Vec<(&'static str, &'static str, &'static str)> {
    vec![
        ("AAA", "K", "a"),
        ("AAC", "N", "b"),
        ("AAG", "K", "c"),
        ("AAT", "N", "d"),
        ("ACA", "T", "e"),
        ("ACC", "T", "f"),
        ("ACG", "T", "g"),
        ("ACT", "T", "h"),
        ("AGA", "R", "i"),
        ("AGC", "S", "j"),
        ("AGG", "R", "k"),
        ("AGT", "S", "l"),
        ("ATA", "I", "m"),
        ("ATC", "I", "n"),
        ("ATG", "M", "o"),
        ("ATT", "I", "p"),
        ("CAA", "Q", "q"),
        ("CAC", "H", "r"),
        ("CAG", "Q", "s"),
        ("CAT", "H", "t"),
        ("CCA", "P", "u"),
        ("CCC", "P", "v"),
        ("CCG", "P", "w"),
        ("CCT", "P", "x"),
        ("CGA", "R", "y"),
        ("CGC", "R", "z"),
        ("CGG", "R", "A"),
        ("CGT", "R", "B"),
        ("CTA", "L", "C"),
        ("CTC", "L", "D"),
        ("CTG", "L", "E"),
        ("CTT", "L", "F"),
        ("GAA", "E", "G"),
        ("GAC", "D", "H"),
        ("GAG", "E", "I"),
        ("GAT", "D", "J"),
        ("GCA", "A", "K"),
        ("GCC", "A", "L"),
        ("GCG", "A", "M"),
        ("GCT", "A", "N"),
        ("GGA", "G", "O"),
        ("GGC", "G", "P"),
        ("GGG", "G", "Q"),
        ("GGT", "G", "R"),
        ("GTA", "V", "S"),
        ("GTC", "V", "T"),
        ("GTG", "V", "U"),
        ("GTT", "V", "V"),
        ("TAA", "*", "W"),
        ("TAC", "Y", "X"),
        ("TAG", "*", "Y"),
        ("TAT", "Y", "Z"),
        ("TCA", "S", "0"),
        ("TCC", "S", "1"),
        ("TCG", "S", "2"),
        ("TCT", "S", "3"),
        ("TGA", "*", "4"),
        ("TGC", "C", "5"),
        ("TGG", "W", "6"),
        ("TGT", "C", "7"),
        ("TTA", "L", "8"),
        ("TTC", "F", "9"),
        ("TTG", "L", "&"),
        ("TTT", "F", "$"),
    ]
}

#[allow(dead_code)]
#[allow(unused_variables)]
pub fn generate_maps(
    translation_table: Vec<(&'static str, &'static str, &'static str)>,
) -> (
    StrStrMap,
    StrStrMap,
    StrStrMap,
    StrStrVecMap,
    StrStrVecMap,
    StrStrVecMap,
) {
    let codon_singular_map = c! {i.0 => i.2, for i in &translation_table};
    let singular_codon_map = c! {i.2 => i.0, for i in &translation_table};
    let codon_amino_map = c! {i.0 => i.1, for i in &translation_table};
    let amino_codon_map = c! {i.1 => c![j.0, for j in &translation_table, if j.1 == i.1],
    for i in &translation_table};
    let amino_singular_map = c! {i.1 => c![j.2, for j in &translation_table, if j.1 == i.1],
    for i in &translation_table};

    let alternatives = c! {i.2 => c![j.2, for j in &translation_table, if j.1 == i.1 && i.2 != j.2],
    for i in &translation_table};

    (
        codon_singular_map,
        singular_codon_map,
        codon_amino_map,
        amino_codon_map,
        amino_singular_map,
        alternatives,
    )
}

use rand::seq::IteratorRandom;
use rand::{thread_rng, Rng};
use std::collections::HashMap;
use std::f64;
use suffix::SuffixTable;

use crate::rayon::iter::IntoParallelRefMutIterator;
use crate::rayon::iter::ParallelIterator;

#[allow(dead_code)]
fn calculate_ars(string: &str, suffix_table: &suffix::SuffixTable) -> f64 {
    // For each position in the string.
    let ars: f64 = string
        .chars()
        .enumerate()
        // Get the per-position score for each substring
        .map(|(idx, _)| {
            let substring = &string[idx..];
            let substring_len = substring.len();
            let mut counter = 0.0;
            for s_idx in 0..substring_len {
                let subsubstring = &substring[..=s_idx];
                if suffix_table.contains(subsubstring) {
                    counter += 1.0;
                } else {
                    break;
                }
            }
            counter
        })
        // Sum them.
        .sum();
    // Do the division and return
    ars / (string.len() as f64)
}

#[derive(Debug, Clone, PartialEq)]
struct Candidate {
    sequence: String,
    fitness: Option<f64>,
}
impl Candidate {
    pub fn new(sequence: String) -> Candidate {
        Candidate {
            sequence,
            fitness: None,
        }
    }

    pub fn mutate(&self, n: usize, alternative_map: &HashMap<&str, Vec<&str>>) -> Candidate {
        // Make a hash set to have chosen positions
        let mut chosen: HashMap<usize, &str> = HashMap::new();
        // Create an RNG instance.
        let mut rng = thread_rng();
        // Choose a number of mutations.
        let no_mutations = rng.gen_range(1, n + 1);
        // Select n positions, and check they have alternatives.
        while chosen.len() != no_mutations {
            // Pick a random position.
            let choice = rng.gen_range(0, self.sequence.len());
            if chosen.contains_key(&choice) {
                continue;
            }
            // Get the corresponding codon encoding.
            let codon_code = self.sequence.chars().nth(choice).unwrap().to_string();
            // Get alternative codon codes.
            let alternative_codes = alternative_map.get::<str>(&codon_code).unwrap();
            // If empty, don't select.
            if alternative_codes.is_empty() {
                continue;
            }
            // Choose an alternative codon, and record that choice.
            chosen.insert(choice, alternative_codes.iter().choose(&mut rng).unwrap());
        }

        // Generate the new sequence
        let new_seq = self
            .sequence
            .chars()
            .enumerate()
            // Use map and match to choose either the new char or the old one.
            .map(|(idx, char)| {
                if chosen.contains_key(&idx) {
                    chosen.get(&idx).unwrap().to_string()
                } else {
                    char.to_string().clone().to_string()
                }
            })
            .collect::<String>();

        // Use the new seq to make a new candidate, and return it.
        Candidate::new(new_seq)
    }

    // Carry out a crossover.
    pub fn crossover(&self, other: &Candidate) -> (Candidate, Candidate) {
        // Create an RNG instance.
        let mut rng = thread_rng();
        // Choose two random positions across the length.
        let mut pos1 = rng.gen_range(0, self.sequence.len());
        let mut pos2 = rng.gen_range(0, self.sequence.len());
        // Ensure that pos1 and pos aren't equal.
        while pos1 == pos2 {
            pos2 = rng.gen_range(0, self.sequence.len());
        }
        // Ensure pos2 < pos1
        if pos1 > pos2 {
            std::mem::swap(&mut pos1, &mut pos2);
        }

        // Create the new sequences.
        let new_seq_one = format!(
            "{}{}{}",
            &self.sequence[..pos1],
            &other.sequence[pos1..pos2],
            &self.sequence[pos2..]
        );
        let new_seq_two = format!(
            "{}{}{}",
            &other.sequence[..pos1],
            &self.sequence[pos1..pos2],
            &other.sequence[pos2..]
        );

        // Generate candidates and return.
        (Candidate::new(new_seq_one), Candidate::new(new_seq_two))
    }

    fn score_min(&mut self, suffix_tables: &[SuffixTable]) {
        let score: f64 = suffix_tables
            .iter()
            .map(|i| calculate_ars(&self.sequence, i))
            .fold(f64::NAN, f64::min);

        self.fitness = Some(score);
    }

    fn score_weighted(&mut self, suffix_tables: &[SuffixTable], weights: &Option<Vec<f64>>) {
        // Calculate ARS for each suffix table, multiply by weight, and sum.
        let mut score = suffix_tables
            .iter()
            .zip(weights.as_ref().unwrap().iter())
            .map(|(st, weight)| weight * calculate_ars(&self.sequence, st))
            .sum();
        // Normalise by the sum of weights.
        score /= weights.as_ref().unwrap().iter().sum::<f64>();
        // Set fitness
        self.fitness = Some(score);
    }

    pub fn score(
        &mut self,
        suffix_tables: &[SuffixTable],
        method: &str,
        weights: &Option<Vec<f64>>,
    ) {
        if let Some(_t) = self.fitness {
            return;
        }

        if method == "weighted" {
            match weights {
                Some(t) => {
                    if t.len() != suffix_tables.len() {
                        panic!("Length of weights and number of organisms given must be equal in weighted mode")
                    }
                }
                None => panic!("Weighted mode requires weights to be given"),
            }
            self.score_weighted(suffix_tables, weights)
        } else if method == "min" {
            self.score_min(suffix_tables)
        } else {
            panic!("Method given for scoring is not supported")
        }
    }
}

fn generate_random_candidates(
    protein: &str,
    map: &HashMap<&str, Vec<&str>>,
    n: usize,
) -> Vec<Candidate> {
    // Generates n random candidates encoding protein using map.
    (0..n)
        .map(|_| {
            // For each character (AA) in the string (protein)
            protein
                .chars()
                .map(|i| {
                    // Create an RNG instance.
                    let mut rng = thread_rng();
                    // Convert character to string
                    let str = i.to_string();
                    // Select a random codon encoding that (AA)
                    let choice = map
                        .get::<str>(&str)
                        .unwrap()
                        .iter()
                        .choose(&mut rng)
                        .unwrap();
                    *choice
                })
                .collect::<String>()
        })
        .map(Candidate::new)
        .collect::<Vec<Candidate>>()
}

fn binary_tournament(population: &mut Vec<Candidate>, target_size: usize) {
    // Make an RNG instance
    let mut rng = thread_rng();
    while population.len() > target_size {
        // Choose two random members of the population.
        let member1_idx = rng.gen_range(0, population.len());
        let mut member2_idx = rng.gen_range(0, population.len());

        while member1_idx == member2_idx {
            member2_idx = rng.gen_range(0, population.len());
        }

        let member1 = population.get(member1_idx).unwrap();
        let member2 = population.get(member2_idx).unwrap();


        let member1_fitness = member1.fitness.unwrap();
        let member2_fitness = member2.fitness.unwrap();

        if member1_fitness > member2_fitness {
            let to_remove = population
                .iter()
                .enumerate()
                .find(|(_, item)| item == &member2)
                .unwrap()
                .0;
            population.remove(to_remove);
        } else {
            let to_remove = population
                .iter()
                .enumerate()
                .find(|(_, item)| item == &member1)
                .unwrap()
                .0;
            population.remove(to_remove);
        }
    }
}

fn crossovers(population: &mut Vec<Candidate>, n_mut: usize) {
    // Make an RNG instance
    let mut rng = thread_rng();
    for _ in 0..n_mut {
        // Choose two random members of the population.
        let member1_idx = rng.gen_range(0, population.len());
        let mut member2_idx = rng.gen_range(0, population.len());

        while member1_idx == member2_idx {
            member2_idx = rng.gen_range(0, population.len());
        }

        let member1 = population.get(member1_idx).unwrap();
        let member2 = population.get(member2_idx).unwrap();

        let (new_a, new_b) = member1.crossover(member2);
        population.push(new_a);
        population.push(new_b)
    }
}

fn mutations(
    population: &mut Vec<Candidate>,
    alternative_map: &HashMap<&str, Vec<&str>>,
    n_cross: usize,
) {
    // Make an RNG instance
    let mut rng = thread_rng();
    for _ in 0..n_cross {
        // Choose a random member of the population.
        let member = population.iter().choose(&mut rng).unwrap();
        let new_a = member.mutate(5, alternative_map);
        population.push(new_a);
    }
}

pub fn run_ea(
    protein: &str,
    suffix_tables: &[SuffixTable],
    amino_singular: &HashMap<&str, Vec<&str>>,
    alternatives: &HashMap<&str, Vec<&str>>,
    n_cross: usize,
    n_mut: usize,
    n_gen: usize,
    n_gen_start: usize,
    method: &str,
    weights: &Option<Vec<f64>>,
) -> (String, f64) {
    // Start by generating random candidates equal to (n_cross * 2) + (n_mut) --> number of new candidates in one generation.
    let mut population = generate_random_candidates(protein, amino_singular, (n_cross * 2) + n_mut);
    population
        .iter_mut()
        .for_each(|i| i.score(suffix_tables, method, weights));

    let highest_fitness = &population
        .iter()
        .map(|i| i.fitness.unwrap())
        .fold(f64::NAN, f64::max);
    info!(
        "Generation 0 : score of best candidate = {}",
        highest_fitness
    );

    // For each of the remaining generations.
    for gen in 1..n_gen {
        // Binary tournament.
        debug!("Generation {}: Carrying out binary tournament", gen);
        binary_tournament(&mut population, n_gen_start);
        debug!("Generation {}: Carrying out crossover events", gen);
        crossovers(&mut population, n_cross);
        debug!("Generation {}: Carrying out mutation events", gen);
        mutations(&mut population, alternatives, n_mut);
        // Score candidates
        debug!("Generation {}: Scoring candidates", gen);
        population
            .par_iter_mut()
            .for_each(|i| i.score(suffix_tables, method, weights));

        // Report
        if gen % 100 == 0 {
            let highest_fitness = &population
                .iter()
                .map(|i| i.fitness.unwrap())
                .fold(f64::NAN, f64::max);
            info!(
                "Generation {} : score of best candidate = {}",
                gen, highest_fitness
            );
        }
    }

    // Get the highest fitness.
    let highest_fitness = &population
        .iter()
        .map(|i| i.fitness.unwrap())
        .fold(f64::NAN, f64::max);
    // Report
    info!(
        "Algorithm complete -- highest fitness {}",
        highest_fitness
    );

    // Get the instance with the highest fitness.
    let best_scoring = c![i, for i in population, if (i.fitness.unwrap() - *highest_fitness).abs() < 1e-12];
    let best = best_scoring.iter().next().unwrap();

    (best.sequence.clone(), best.fitness.unwrap())
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    use suffix::SuffixTable;

    #[test]
    fn test_ars_v1() {
        let st = SuffixTable::new("ACTG");
        let string = String::from("ACTG");
        assert_eq!(calculate_ars(&string, &st), 2.5);
    }
    #[test]
    fn test_ars_v2() {
        let st = SuffixTable::new("ACTG");
        let string = String::from("ACTA");
        assert_eq!(calculate_ars(&string, &st), 1.75);
    }
    #[test]
    fn test_ars_v3() {
        let st = SuffixTable::new("ACTG");
        let string = String::from("AAAA");
        assert_eq!(calculate_ars(&string, &st), 1.0);
    }
    #[test]
    fn test_ars_v4() {
        let st = SuffixTable::new("ACTG");
        let string = String::from("ACAC");
        assert_eq!(calculate_ars(&string, &st), 1.5);
    }
}

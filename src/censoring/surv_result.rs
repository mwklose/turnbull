use std::iter::zip;

use nalgebra::DVector;

pub struct SurvResult {
    pub times: Vec<(f64, f64, f64)>,
}

impl SurvResult {
    pub fn new_result(times_vec: Vec<f64>, density_vec: DVector<f64>) -> Self {
        let mut result_vec = Vec::new();

        if times_vec.len() != density_vec.len() {
            panic!("Length of unique event times does not equal density vector.");
        }

        let mut cumsum = density_vec.clone();
        cumsum.iter_mut().fold(0.0, |acc, x| {
            *x += acc;
            *x
        });

        zip(times_vec.iter(), cumsum.iter()).for_each(|(&t, &d)| {
            result_vec.push((t, 1.0 - d, 0.0));
        });

        SurvResult { times: result_vec }
    }

    pub fn get_survival_table(&self) -> Vec<(f64, f64, f64)> {
        self.times.clone()
    }
}

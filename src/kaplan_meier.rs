use core::f64;
use std::ops::Div;

use nalgebra::{DMatrix, DVector, Dyn, Matrix, VecStorage, Vector, U1};

use crate::censoring::{surv_result::SurvResult, Surv};
extern crate nalgebra as na;

pub fn kaplan_meier(events: Vec<Surv>, weights: Option<Vec<f64>>) -> SurvResult {
    let mut unique_event_times = events
        .iter()
        .filter(|&ev| !ev.is_censored() && !ev.is_interval_censored())
        .map(|ev| ev.get_exit_time())
        .collect::<Vec<f64>>();

    unique_event_times.dedup();
    unique_event_times.sort_by(|a, b| a.total_cmp(b));
    unique_event_times.push(f64::INFINITY);

    let mut density_vec = DVector::from_vec(unique_event_times.clone())
        .map(|_| 1.0 / unique_event_times.len() as f64);

    let weight_matrix = match weights {
        Some(x) => DVector::from_vec(x),
        None => DVector::from_element(events.len(), 1_f64),
    };

    // TODO: get indicators of censoring set at current time point
    let dm = DMatrix::from_fn(events.len(), unique_event_times.len(), |i, j| {
        let surv = events.get(i).unwrap();
        let event_time = unique_event_times.get(j).unwrap();
        if surv.get_exit_l() <= *event_time && surv.get_exit_time() >= *event_time {
            return 1.0;
        }

        return 0.0;
    });

    // TODO: loop here
    let mut tol = 1e9;
    while tol > 1e-9 {
        (tol, density_vec) = km_helper(&dm, &density_vec, &weight_matrix);
    }

    // TODO: variance estimation

    return SurvResult::new_result(unique_event_times, density_vec);
}

fn km_helper(
    event_matrix_set: &Matrix<f64, Dyn, Dyn, VecStorage<f64, Dyn, Dyn>>,
    density_vec: &Matrix<f64, Dyn, U1, VecStorage<f64, Dyn, U1>>,
    weights_matrix: &Vector<f64, Dyn, VecStorage<f64, Dyn, U1>>,
) -> (f64, Matrix<f64, Dyn, U1, VecStorage<f64, Dyn, U1>>) {
    // Expand from px1 matrix into nxp matrix
    let nrows = event_matrix_set.nrows();
    let intermediates = (density_vec * DMatrix::from_element(1, nrows, 1.0)).transpose();

    // Expectation: allocate densities from censored individuals
    // So then can do component-wise multiplication
    //
    let mut component_wise = event_matrix_set.component_mul(&intermediates);

    // Noramlize rows so they sum up to one
    for mut row in component_wise.row_iter_mut() {
        let row_sum = row.sum();
        row /= row_sum
    }

    // Get column sums

    let mut colsums = component_wise.tr_mul(weights_matrix);

    println!("Colsums: {:?}\n", colsums);
    // Maximization
    // Normalize column sums

    let total_sum = colsums.sum();
    colsums.iter_mut().for_each(|c| {
        *c = c.div(total_sum);
    });

    println!("Result: {:?}\n", colsums);
    (colsums.metric_distance(&density_vec), colsums)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn km_simple() {
        let surv_vec = vec![
            Surv::new_event(1, 0.0, 2.0),
            Surv::new_event(1, 0.0, 4.0),
            Surv::new_event(1, 0.0, 6.0),
            Surv::new_event(1, 0.0, 8.0),
            Surv::new_event(1, 0.0, 10.0),
        ];

        let result = kaplan_meier(surv_vec, None);

        let expected = vec![
            (2.0, 0.8, 0.0),
            (4.0, 0.6, 0.0),
            (6.0, 0.4, 0.0),
            (8.0, 0.2, 0.0),
            (10.0, 0.0, 0.0),
            (f64::INFINITY, 0.0, 0.0),
        ];

        for (i, t) in result.get_survival_table().iter().enumerate() {
            println!("Time {}: {:?}", i, t);

            assert!(t.0 == expected.get(i).unwrap().0);
            assert!((t.1 - expected.get(i).unwrap().1) < 1e-9);
            assert!(t.2 == expected.get(i).unwrap().2);
        }
    }

    #[test]
    fn km_simple_censoring() {
        let surv_vec = vec![
            Surv::new_event(1, 0.0, 2.0),
            Surv::new_event(1, 0.0, 4.0),
            Surv::new_censor(0.0, 6.0),
            Surv::new_event(1, 0.0, 8.0),
            Surv::new_event(1, 0.0, 10.0),
        ];

        let result = kaplan_meier(surv_vec, None);

        let expected = vec![
            (2.0, 0.8, 0.0),
            (4.0, 0.6, 0.0),
            (8.0, 0.3, 0.0),
            (10.0, 0.0, 0.0),
            (f64::INFINITY, 0.0, 0.0),
        ];

        for (i, t) in result.get_survival_table().iter().enumerate() {
            println!("Time {}: {:?}", i, t);

            assert!(t.0 == expected.get(i).unwrap().0);
            assert!((t.1 - expected.get(i).unwrap().1) < 1e-9);
            assert!(t.2 == expected.get(i).unwrap().2);
        }
    }
}

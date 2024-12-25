use core::f64;

pub mod surv_result;

#[derive(Debug)]
pub struct Surv {
    indicator: usize,
    entry_time: f64,
    exit_l: f64,
    exit_r: f64,
}

impl Surv {
    pub fn new(ind: usize, entry_time: f64, exit_l: f64, exit_r: f64) -> Surv {
        Surv {
            indicator: ind,
            entry_time,
            exit_l,
            exit_r,
        }
    }

    pub fn new_event(ind: usize, entry_time: f64, event_time: f64) -> Surv {
        Surv {
            indicator: ind,
            entry_time,
            exit_l: event_time,
            exit_r: event_time,
        }
    }

    pub fn new_censor(entry_time: f64, event_time: f64) -> Surv {
        Surv {
            indicator: 0,
            entry_time,
            exit_l: event_time,
            exit_r: f64::INFINITY,
        }
    }

    pub fn from_exit_vec(vec: Vec<(f64, usize)>) -> Vec<Surv> {
        vec.iter()
            .map(|(event_time, event_ind)| Surv {
                indicator: *event_ind,
                entry_time: 0.0,
                exit_l: *event_time,
                exit_r: *event_time,
            })
            .collect::<Vec<Surv>>()
    }

    pub fn is_censored(&self) -> bool {
        self.indicator == 0
    }

    pub fn is_interval_censored(&self) -> bool {
        !(self.exit_l == self.exit_r)
    }

    pub fn get_entry_time(&self) -> f64 {
        self.entry_time
    }

    pub fn get_exit_l(&self) -> f64 {
        self.exit_l
    }

    pub fn get_exit_time(&self) -> f64 {
        self.exit_r
    }
}

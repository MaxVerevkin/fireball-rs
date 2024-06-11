use std::fmt::Write;
use std::io;
use std::{borrow::Borrow, collections::HashMap};

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use common::{
    structs::{Line, UnitVec3, Vec3},
    Sigmas,
};

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct Db {
    pub events: HashMap<String, Event>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Event {
    pub answer: Answer,
    pub observers_count: u32,
    pub avg_observer: Vec3,
    pub runs: HashMap<ParamsKey, Vec<Stage>>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct Answer {
    pub point: Vec3,
    pub dir: Option<UnitVec3>,
}

#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[serde(transparent)]
pub struct ParamsKey(pub String);

impl Borrow<str> for ParamsKey {
    fn borrow(&self) -> &str {
        &self.0
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Stage {
    pub title: String,
    pub extra_info: Option<String>,
    pub traj: Line,
    pub sigmas: Option<Sigmas>,
}

impl Db {
    pub fn read(path: Option<&str>) -> Result<Self> {
        match std::fs::File::open(path.unwrap_or("db.json")) {
            Ok(file) => {
                let reader = io::BufReader::new(file);
                let db = serde_json::from_reader(reader).context("could not deserialize DB")?;
                Ok(db)
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(Self::default()),
            Err(e) => Err(anyhow::Error::from(e).context("could not open db.json")),
        }
    }

    pub fn write(&self) -> Result<()> {
        let file = std::fs::File::create("db.json").context("could not open db.json")?;
        let writer = io::BufWriter::new(file);
        serde_json::to_writer(writer, self).context("could not serialize DB")?;
        Ok(())
    }

    pub fn add_run(
        &mut self,
        data: &common::obs_data::Data,
        params: crate::solver::Params,
        stages: Vec<Stage>,
    ) {
        let event_name = data.name.clone().unwrap();
        let answer = data.answer.unwrap();
        let observers_count = data.samples.len() as u32;

        let answer = Answer {
            point: answer.point.unwrap(),
            dir: answer.traj.map(|traj| traj.direction),
        };
        let entry = self.events.entry(event_name).or_insert_with(|| Event {
            answer,
            observers_count,
            avg_observer: data
                .samples
                .iter()
                .map(|s| s.location)
                .fold(Vec3::new(0.0, 0.0, 0.0), |a, b| a + b)
                / data.samples.len() as f64,
            runs: HashMap::new(),
        });
        entry.answer = answer;
        entry.observers_count = observers_count;
        entry.runs.insert(ParamsKey::from_params(params), stages);
    }
}

impl ParamsKey {
    fn from_params(params: crate::solver::Params) -> Self {
        let mut res = String::new();
        write!(res, "da_k={};", params.da_k).unwrap();
        write!(res, "az_k={};", params.az_k).unwrap();
        if params.correct_altitudes {
            res.push_str("correct_altitudes;");
        }
        if params.no_altitudes {
            res.push_str("no_altitudes;");
        }
        if params.no_azimuths {
            res.push_str("no_azimuths;");
        }
        if params.no_da_flip {
            res.push_str("no_da_flip;");
        }
        Self(res)
    }
}

use std::collections::HashMap;
use std::fmt::Write;
use std::io;

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
    pub runs: HashMap<ParamsKey, Vec<Stage>>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Answer {
    pub point: Vec3,
    pub dir: Option<UnitVec3>,
}

#[derive(Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
#[serde(transparent)]
pub struct ParamsKey(pub String);

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Stage {
    pub title: String,
    pub extra_info: Option<String>,
    pub traj: Line,
    pub sigmas: Option<Sigmas>,
}

impl Db {
    pub fn read() -> Result<Self> {
        match std::fs::File::open("db.json") {
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
        event_name: String,
        answer: common::obs_data::Answer,
        params: crate::solver::Params,
        stages: Vec<Stage>,
    ) {
        let answer = Answer {
            point: answer.point.unwrap(),
            dir: answer.traj.map(|traj| traj.direction),
        };
        self.events
            .entry(event_name)
            .or_insert_with(|| Event {
                answer,
                runs: HashMap::new(),
            })
            .runs
            .insert(ParamsKey::from_params(params), stages);
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
        Self(res)
    }
}

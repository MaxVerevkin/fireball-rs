#!/bin/sh
date >> log.txt
cargo r --release -- "$@" | tee -a log.txt

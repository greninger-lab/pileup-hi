// Copyright 2021 Olivier Kannengieser
//
// Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

#![cfg_attr(feature="nightly",feature(core_intrinsics))]

use likely_stable::{likely,unlikely,if_likely,if_unlikely};

#[test]
fn test_likely() {
    assert!(if likely(true) { true} else {false});
    assert!(if likely(false) { false} else {true});
}
#[test]
fn test_unlikely() {
    assert!(if unlikely(true) { true} else {false});
    assert!(if unlikely(false) { false} else {true});
}
#[test]
fn test_if_likely () {

    let mut taken = false;
    if_likely!(0<1 => {taken=true});
    assert!(taken);

    let mut taken = false;
    if_likely!(1<0 => {taken=true});
    assert!(!taken);

    assert!(if_likely!(0<1 => {true} else {false}));

    assert!(if_likely!(0>1 => {false} else {true}));

    let mut taken = false;
    if_likely!(let Some(_v) = Some(0) => {taken=true});
    assert!(taken);

    let mut taken = false;
    if_likely!(let Some(_v) = Option::<i32>::None => {taken=true});
    assert!(!taken);

    assert!(if_likely!(let Some(_v) = Some(0) => {true} else {false}));

    assert!(if_likely!(let Some(_v) = Option::<i32>::None => {false} else {true}));
}
#[test]
fn test_if_unlikely () {

    let mut taken = false;
    if_unlikely!(0<1 => {taken=true});
    assert!(taken);

    let mut taken = false;
    if_unlikely!(1<0 => {taken=true});
    assert!(!taken);

    assert!(if_unlikely!(0<1 => {true} else {false}));

    assert!(if_unlikely!(0>1 => {false} else {true}));

    let mut taken = false;
    if_unlikely!(let Some(_v) = Some(0) => {taken=true});
    assert!(taken);

    let mut taken = false;
    if_unlikely!(let Some(_v) = Option::<i32>::None => {taken=true});
    assert!(!taken);

    assert!(if_unlikely!(let Some(_v) = Some(0) => {true} else {false}));

    assert!(if_unlikely!(let Some(_v) = Option::<i32>::None => {false} else {true}));
}



// Copyright 2021 Olivier Kannengieser
//
// Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
// http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
// http://opensource.org/licenses/MIT>, at your option. This file may not be
// copied, modified, or distributed except according to those terms.

#![no_std]
#![cfg_attr(nightly, feature(core_intrinsics))]

//! This crates brings [likely](core::intrinsics::likely) 
//! and [unlikely](core::intrinsics::unlikely) branch prediction hints to stable rust
//!
//! Likelyness hint only affects static branch prediction. So a good likelyness indication
//! is the one that indicates the likelyness a branch will be taken when the condition has not been executed
//! since a relatively long time. Profile measurments are very bad indicator on weither to use likely hint
//! or unlikely hint for a given branching statement.
//!
//! The fact it only affects the static branch predictor makes likely-hood optimizations extremely difficult to benchmark.
//! Indeed, all micro-benchmark needs to loop over the benchmarked code in order to
//! have an accurate timing. But since we loop, the branch predictor use its dynamic prediction
//! ability and the likely-hood ness hint is not used. 
//! ```
//! use likely_stable::{likely,unlikely};
//! use rand::random;
//!
//! if likely(random::<i32>() > 10) {
//!     println!("likely!")
//! } else {
//!     println!("unlikely!")
//! }
//! ```
//!
//! It also provides [macro@if_likely] and [macro@if_unlikely] for branch prediction
//! for `if let` statements.
//! ```
//! use likely_stable::if_likely;
//! use rand::random;
//!
//! let v = Some(random()).filter(|v:&i32| *v > 10);
//!
//! if_likely!{let Some(v) = v => {
//!     println!("likely!")
//! } else {
//!     println!("unlikely!")
//! }};
//! ```
//!
//! Moreover traits [LikelyBool], [LikelyOption] and [LikelyResult] provides *likely*
//! and *unlikely* versions of the methods commonly used for types [bool], [Option] and
//! [Result]
//! ```
//! use likely_stable::LikelyOption;
//! use rand::random;
//!
//! let v = Some(random()).filter(|v:&i32| *v > 10);
//!
//! v.map_or_else_likely(
//!     || println!("unlikely"),
//!     |v| println!("likely {}",v));
//! ```

#[cfg(feature = "nightly")]
pub use core::intrinsics::{likely, unlikely};

use core::hint::unreachable_unchecked;

#[cfg(not(feature = "nightly"))]
#[inline(always)]
/// Brings [likely](core::intrinsics::likely) to stable rust.
pub const fn likely(b: bool) -> bool {
    #[allow(clippy::needless_bool)]
    if (1i32).checked_div(if b { 1 } else { 0 }).is_some() {
        true
    } else {
        false
    }
}

#[cfg(not(feature = "nightly"))]
#[inline(always)]
/// Brings [unlikely](core::intrinsics::unlikely) to stable rust.
pub const fn unlikely(b: bool) -> bool {
    #[allow(clippy::needless_bool)]
    if (1i32).checked_div(if b { 0 } else { 1 }).is_none() {
        true
    } else {
        false
    }
}

/// Likely trait for bool
///
/// `likely` method suffix means *likely true*.
pub trait LikelyBool: Sized {
    fn then_likely<T, F: FnOnce() -> T>(self, f: F) -> Option<T>;
    fn then_unlikely<T, F: FnOnce() -> T>(self, f: F) -> Option<T>;
}

/// Likely trait for Options
///
/// `likely` method suffix means *likely Some*.
pub trait LikelyOption: Sized {
    type Value;
    fn and_likely<U>(self, res: Option<U>) -> Option<U>;
    fn and_unlikely<U>(self, res: Option<U>) -> Option<U>;
    fn and_then_likely<U, F: FnOnce(Self::Value) -> Option<U>>(self, f: F) -> Option<U>;
    fn and_then_unlikely<U, F: FnOnce(Self::Value) -> Option<U>>(self, f: F) -> Option<U>;
    fn get_or_insert_likely(&mut self, value: Self::Value) -> &mut Self::Value;
    fn get_or_insert_unlikely(&mut self, value: Self::Value) -> &mut Self::Value;
    fn get_or_insert_with_likely<F: FnOnce() -> Self::Value>(&mut self, f: F) -> &mut Self::Value;
    fn get_or_insert_with_unlikely<F: FnOnce() -> Self::Value>(&mut self, f: F)
        -> &mut Self::Value;
    fn filter_likely<P: FnOnce(&Self::Value) -> bool>(self, predicate: P) -> Self;
    fn filter_unlikely<P: FnOnce(&Self::Value) -> bool>(self, predicate: P) -> Self;
    fn map_likely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Option<U>;
    fn map_unlikely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Option<U>;
    fn map_or_likely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U;
    fn map_or_unlikely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U;
    fn map_or_else_likely<U, D: FnOnce() -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U;
    fn map_or_else_unlikely<U, D: FnOnce() -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U;
    fn or_likely(self, b: Self) -> Self;
    fn or_unlikely(self, b: Self) -> Self;
    fn or_else_likely<F: FnOnce() -> Self>(self, b: F) -> Self;
    fn or_else_unlikely<F: FnOnce() -> Self>(self, b: F) -> Self;
    fn unwrap_or_likely(self, d: Self::Value) -> Self::Value;
    fn unwrap_or_unlikely(self, d: Self::Value) -> Self::Value;
    fn unwrap_or_else_likely<F: FnOnce() -> Self::Value>(self, b: F) -> Self::Value;
    fn unwrap_or_else_unlikely<F: FnOnce() -> Self::Value>(self, b: F) -> Self::Value;
}
/// Likely trait for Result
///
/// `likely` method suffix means *likely Ok*.
pub trait LikelyResult: Sized {
    type Value;
    type Error;
    fn and_likely<U>(self, res: Result<U, Self::Error>) -> Result<U, Self::Error>;
    fn and_unlikely<U>(self, res: Result<U, Self::Error>) -> Result<U, Self::Error>;
    fn and_then_likely<U, F: FnOnce(Self::Value) -> Result<U, Self::Error>>(
        self,
        f: F,
    ) -> Result<U, Self::Error>;
    fn and_then_unlikely<U, F: FnOnce(Self::Value) -> Result<U, Self::Error>>(
        self,
        f: F,
    ) -> Result<U, Self::Error>;
    fn map_likely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Result<U, Self::Error>;
    fn map_unlikely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Result<U, Self::Error>;
    fn map_err_likely<U, F: FnOnce(Self::Error) -> U>(self, f: F) -> Result<Self::Value, U>;
    fn map_err_unlikely<U, F: FnOnce(Self::Error) -> U>(self, f: F) -> Result<Self::Value, U>;
    fn map_or_likely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U;
    fn map_or_unlikely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U;
    fn map_or_else_likely<U, D: FnOnce(Self::Error) -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U;
    fn map_or_else_unlikely<U, D: FnOnce(Self::Error) -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U;
    fn or_likely<F>(self, res: Result<Self::Value, F>) -> Result<Self::Value, F>;
    fn or_unlikely<F>(self, res: Result<Self::Value, F>) -> Result<Self::Value, F>;
    fn or_else_likely<F: FnOnce(Self::Error) -> Result<Self::Value, F>>(
        self,
        op: F,
    ) -> Result<Self::Value, F>;
    fn or_else_unlikely<F: FnOnce(Self::Error) -> Result<Self::Value, F>>(
        self,
        op: F,
    ) -> Result<Self::Value, F>;
    fn unwrap_or_likely(self, default: Self::Value) -> Self::Value;
    fn unwrap_or_unlikely(self, default: Self::Value) -> Self::Value;
    fn unwrap_or_else_likely<F: FnOnce(Self::Error) -> Self::Value>(self, op: F) -> Self::Value;
    fn unwrap_or_else_unlikely<F: FnOnce(Self::Error) -> Self::Value>(self, op: F) -> Self::Value;
}

/// If statement which inform the compiler
/// that the first branch will be the most taken branch
///
/// It is usefull as a `if let` statment:
/// ```
/// use likely_stable::if_likely;
/// use rand::random;
/// 
/// let v = Some(random()).filter(|v:&i32| *v > 10);
/// 
/// if_likely!{let Some(v) = v => {
///     println!("likely!")
/// } else {
///     println!("unlikely!")
/// }};
/// ```
#[macro_export]
macro_rules! if_likely {
    ($cond:expr => $e:block) => {
       {
       if $crate::likely($cond)
               $e
       }
    };
    ($cond:expr => $e:block else $o:block) => {
       {
       if $crate::likely($cond)
           $e
        else
           $o
        }
    };
    (let $v:pat = $cond:expr => $e:block) => {
       {
       let __cond_expr = $cond;
       if $crate::likely(if let $v = &__cond_expr { true } else {false}) {
           if let $v = __cond_expr
               $e
            else {
               //SAFETY: pattern matching should not mutate __cond_expr
               //  so the too if let shall match equaly
               unsafe{::core::hint::unreachable_unchecked()}
           }
       }
       }
    };
    (let $v:pat = $cond:expr => $e:block else $o:block) => {
       {
       let __cond_expr = $cond;
       if $crate::likely( if  let $v = &__cond_expr { true } else {false}) {
           if let $v = __cond_expr
               $e
            else {
               //SAFETY: pattern matching should not mutate __cond_expr
               //  so the too if let shall match equaly
               unsafe{::core::hint::unreachable_unchecked()}
           }
       } else
           $o
        }
    }
}

/// If statement which inform the compiler
/// that the first branch will be the most taken branch
///
/// It is usefull as a `if let` statment:
/// ```
/// use likely_stable::if_unlikely;
/// use rand::random;
/// 
/// let v = Some(random()).filter(|v:&i32| *v > 10);
/// 
/// if_unlikely!{let None = v => {
///     println!("unlikely!")
/// } else {
///     println!("likely!")
/// }};
/// ```
#[macro_export]
macro_rules! if_unlikely {
    ($cond:expr => $e:block) => {
       {
       if $crate::unlikely($cond)
           $e
       }
    };
    ($cond:expr => $e:block else $o:block) => {
       {
       if $crate::unlikely($cond)
           $e
        else
           $o
        }
    };
    (let $v:pat = $cond:expr => $e:block) => {
       {
       let __cond_expr = $cond;
       if $crate::unlikely(if let $v = &__cond_expr { true } else {false}) {
           if let $v = __cond_expr
               $e
            else {
               //SAFETY: pattern matching should not mutate __cond_expr
               //  so the too if let shall match equaly
               unsafe{::core::hint::unreachable_unchecked()}
           }
       }
       }
    };
    (let $v:pat = $cond:expr => $e:block else $o:block) => {
       {
       let __cond_expr = $cond;
       if $crate::unlikely( if  let $v = &__cond_expr { true } else {false}) {
           if let $v = __cond_expr
               $e
            else {
               //SAFETY: pattern matching should not mutate __cond_expr
               //  so the too if let shall match equaly
               unsafe{::core::hint::unreachable_unchecked()}
           }
       } else
           $o
        }
    }
}

impl LikelyBool for bool {
    #[inline]
    fn then_likely<T, F: FnOnce() -> T>(self, f: F) -> Option<T> {
        if likely(self) {
            Some(f())
        } else {
            None
        }
    }

    #[inline]
    fn then_unlikely<T, F: FnOnce() -> T>(self, f: F) -> Option<T> {
        if unlikely(self) {
            Some(f())
        } else {
            None
        }
    }
}

#[inline(always)]
fn apply_opt_likely<V, R, FOk: FnOnce(V) -> R, FErr: FnOnce() -> R>(
    this: Option<V>,
    fok: FOk,
    ferr: FErr,
) -> R {
    if likely(this.is_some()) {
        //SAFETY: is_some() does not mutate the option 
        match this {
            Some(v) => fok(v),
            None => unsafe { unreachable_unchecked() },
        }
    } else {
        //SAFETY: is_some() does not mutate the option 
        match this {
            None => ferr(),
            Some(_) => unsafe { unreachable_unchecked() },
        }
    }
}
#[inline(always)]
fn apply_opt_unlikely<V, R, FOk: FnOnce(V) -> R, FErr: FnOnce() -> R>(
    this: Option<V>,
    fok: FOk,
    ferr: FErr,
) -> R {
    if unlikely(this.is_some()) {
        //SAFETY: is_some() does not mutate the option 
        match this {
            Some(v) => fok(v),
            None => unsafe { unreachable_unchecked() },
        }
    } else {
        //SAFETY: is_some() does not mutate the option 
        match this {
            None => ferr(),
            Some(_) => unsafe { unreachable_unchecked() },
        }
    }
}

impl<T> LikelyOption for Option<T> {
    type Value = T;

    #[inline(always)]
    fn and_likely<U>(self, res: Option<U>) -> Option<U> {
        apply_opt_likely(self, |_| res, || None)
    }
    #[inline(always)]
    fn and_unlikely<U>(self, res: Option<U>) -> Option<U> {
        apply_opt_unlikely(self, |_| res, || None)
    }
    #[inline]
    fn and_then_likely<U, F: FnOnce(Self::Value) -> Option<U>>(self, f: F) -> Option<U> {
        apply_opt_likely(self, f, || None)
    }
    #[inline]
    fn and_then_unlikely<U, F: FnOnce(Self::Value) -> Option<U>>(self, f: F) -> Option<U> {
        apply_opt_unlikely(self, f, || None)
    }
    #[inline]
    fn filter_likely<P: FnOnce(&T) -> bool>(self, predicate: P) -> Option<T> {
        apply_opt_likely(
            self,
            |v| if predicate(&v) { Some(v) } else { None },
            || None,
        )
    }
    #[inline]
    fn filter_unlikely<P: FnOnce(&T) -> bool>(self, predicate: P) -> Option<T> {
        apply_opt_unlikely(
            self,
            |v| if predicate(&v) { Some(v) } else { None },
            || None,
        )
    }
    #[inline(always)]
    fn get_or_insert_likely(&mut self, value: Self::Value) -> &mut Self::Value {
        self.get_or_insert_with_likely(|| value)
    }
    #[inline(always)]
    fn get_or_insert_unlikely(&mut self, value: Self::Value) -> &mut Self::Value {
        self.get_or_insert_with_unlikely(|| value)
    }
    #[inline]
    fn get_or_insert_with_likely<F: FnOnce() -> T>(&mut self, f: F) -> &mut Self::Value {
        if likely(self.is_some()) {
            //SAFETY: is_some() does not mutate the option 
            match self {
                Some(v) => v,
                None => unsafe { unreachable_unchecked() },
            }
        } else {
            match self {
                None => {
                    *self = Some(f());
                    //SAFETY: is_some() does not mutate the option 
                    match self {
                        Some(v) => v,
                        None => unsafe { unreachable_unchecked() },
                    }
                }
                //SAFETY: is_some() does not mutate the option 
                Some(_) => unsafe { unreachable_unchecked() },
            }
        }
    }
    #[inline]
    fn get_or_insert_with_unlikely<F: FnOnce() -> T>(&mut self, f: F) -> &mut Self::Value {
        if unlikely(self.is_some()) {
            //SAFETY: is_some() does not mutate the option 
            match self {
                Some(v) => v,
                None => unsafe { unreachable_unchecked() },
            }
        } else {
            //SAFETY: is_some() does not mutate the option 
            match self {
                None => {
                    *self = Some(f());
                    match self {
                        Some(v) => v,
                        None => unsafe { unreachable_unchecked() },
                    }
                }
                Some(_) => unsafe { unreachable_unchecked() },
            }
        }
    }
    #[inline]
    fn map_likely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Option<U> {
        apply_opt_likely(self, |v| Some(f(v)), || None)
    }

    #[inline]
    fn map_unlikely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Option<U> {
        apply_opt_unlikely(self, |v| Some(f(v)), || None)
    }
    #[inline]
    fn map_or_likely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U {
        apply_opt_likely(self, f, || default)
    }
    #[inline]
    fn map_or_unlikely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U {
        apply_opt_unlikely(self, f, || default)
    }
    #[inline]
    fn map_or_else_likely<U, D: FnOnce() -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U {
        apply_opt_likely(self, f, default)
    }
    #[inline]
    fn map_or_else_unlikely<U, D: FnOnce() -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U {
        apply_opt_unlikely(self, f, default)
    }
    #[inline(always)]
    fn or_likely(self, b: Self) -> Self {
        apply_opt_likely(self, |v| Some(v), || b)
    }
    #[inline(always)]
    fn or_unlikely(self, b: Self) -> Self {
        apply_opt_unlikely(self, |v| Some(v), || b)
    }
    #[inline]
    fn or_else_likely<F: FnOnce() -> Self>(self, b: F) -> Self {
        apply_opt_likely(self, |v| Some(v), b)
    }
    #[inline]
    fn or_else_unlikely<F: FnOnce() -> Self>(self, b: F) -> Self {
        apply_opt_unlikely(self, |v| Some(v), b)
    }
    #[inline(always)]
    fn unwrap_or_likely(self, d: Self::Value) -> Self::Value {
        apply_opt_likely(self, |v| v, || d)
    }
    #[inline(always)]
    fn unwrap_or_unlikely(self, d: Self::Value) -> Self::Value {
        apply_opt_unlikely(self, |v| v, || d)
    }
    #[inline]
    fn unwrap_or_else_likely<F: FnOnce() -> Self::Value>(self, b: F) -> Self::Value {
        apply_opt_likely(self, |v| v, b)
    }
    #[inline]
    fn unwrap_or_else_unlikely<F: FnOnce() -> Self::Value>(self, b: F) -> Self::Value {
        apply_opt_unlikely(self, |v| v, b)
    }
}

#[inline(always)]
fn apply_res_likely<V, E, R, FOk: FnOnce(V) -> R, FErr: FnOnce(E) -> R>(
    this: Result<V, E>,
    fok: FOk,
    ferr: FErr,
) -> R {
    if likely(this.is_ok()) {
        //SAFETY: is_ok() does not mutate the option 
        match this {
            Ok(v) => fok(v),
            Err(_) => unsafe { unreachable_unchecked() },
        }
    } else {
        //SAFETY: is_ok() does not mutate the option 
        match this {
            Err(e) => ferr(e),
            Ok(_) => unsafe { unreachable_unchecked() },
        }
    }
}
#[inline(always)]
fn apply_res_unlikely<V, E, R, FOk: FnOnce(V) -> R, FErr: FnOnce(E) -> R>(
    this: Result<V, E>,
    fok: FOk,
    ferr: FErr,
) -> R {
    if unlikely(this.is_ok()) {
        //SAFETY: is_ok() does not mutate the option 
        match this {
            Ok(v) => fok(v),
            Err(_) => unsafe { unreachable_unchecked() },
        }
    } else {
        //SAFETY: is_ok() does not mutate the option 
        match this {
            Err(e) => ferr(e),
            Ok(_) => unsafe { unreachable_unchecked() },
        }
    }
}

impl<V, E> LikelyResult for Result<V, E> {
    type Value = V;
    type Error = E;

    #[inline(always)]
    fn and_likely<U>(self, res: Result<U, E>) -> Result<U, E> {
        apply_res_likely(self, |_| res, |e| Err(e))
    }
    #[inline(always)]
    fn and_unlikely<U>(self, res: Result<U, E>) -> Result<U, E> {
        apply_res_unlikely(self, |_| res, |e| Err(e))
    }
    #[inline]
    fn and_then_likely<U, F: FnOnce(Self::Value) -> Result<U, E>>(self, f: F) -> Result<U, E> {
        apply_res_likely(self, f, |e| Err(e))
    }
    #[inline]
    fn and_then_unlikely<U, F: FnOnce(Self::Value) -> Result<U, E>>(self, f: F) -> Result<U, E> {
        apply_res_unlikely(self, f, |e| Err(e))
    }
    #[inline]
    fn map_likely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Result<U, E> {
        apply_res_likely(self, |v| Ok(f(v)), |e| Err(e))
    }
    #[inline]
    fn map_unlikely<U, F: FnOnce(Self::Value) -> U>(self, f: F) -> Result<U, E> {
        apply_res_unlikely(self, |v| Ok(f(v)), |e| Err(e))
    }
    #[inline]
    fn map_err_likely<U, F: FnOnce(E) -> U>(self, f: F) -> Result<V, U> {
        apply_res_likely(self, |v| Ok(v), |e| Err(f(e)))
    }
    #[inline]
    fn map_err_unlikely<U, F: FnOnce(E) -> U>(self, f: F) -> Result<V, U> {
        apply_res_unlikely(self, |v| Ok(v), |e| Err(f(e)))
    }
    #[inline]
    fn map_or_likely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U {
        apply_res_likely(self, f, |_| default)
    }
    #[inline]
    fn map_or_unlikely<U, F: FnOnce(Self::Value) -> U>(self, default: U, f: F) -> U {
        apply_res_unlikely(self, f, |_| default)
    }
    #[inline]
    fn map_or_else_likely<U, D: FnOnce(Self::Error) -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U {
        apply_res_likely(self, f, default)
    }
    #[inline]
    fn map_or_else_unlikely<U, D: FnOnce(Self::Error) -> U, F: FnOnce(Self::Value) -> U>(
        self,
        default: D,
        f: F,
    ) -> U {
        apply_res_unlikely(self, f, default)
    }
    #[inline(always)]
    fn unwrap_or_likely(self, default: V) -> V {
        apply_res_likely(self, |v| v, |_| default)
    }
    #[inline(always)]
    fn unwrap_or_unlikely(self, default: V) -> V {
        apply_res_unlikely(self, |v| v, |_| default)
    }
    #[inline]
    fn unwrap_or_else_likely<F: FnOnce(Self::Error) -> Self::Value>(self, op: F) -> V {
        apply_res_likely(self, |v| v, op)
    }
    #[inline]
    fn unwrap_or_else_unlikely<F: FnOnce(Self::Error) -> Self::Value>(self, op: F) -> V {
        apply_res_unlikely(self, |v| v, op)
    }
    #[inline]
    fn or_likely<F>(self, res: Result<Self::Value, F>) -> Result<Self::Value, F> {
        apply_res_likely(self, |v| Ok(v), |_| res)
    }
    #[inline]
    fn or_unlikely<F>(self, res: Result<Self::Value, F>) -> Result<Self::Value, F> {
        apply_res_unlikely(self, |v| Ok(v), |_| res)
    }
    #[inline]
    fn or_else_likely<F: FnOnce(Self::Error) -> Result<Self::Value, F>>(
        self,
        op: F,
    ) -> Result<Self::Value, F> {
        apply_res_likely(self, |v| Ok(v), op)
    }
    #[inline]
    fn or_else_unlikely<F: FnOnce(Self::Error) -> Result<Self::Value, F>>(
        self,
        op: F,
    ) -> Result<Self::Value, F> {
        apply_res_unlikely(self, |v| Ok(v), op)
    }
}

#[cfg(feature = "check_assembly")]
/// # How to check the resulting assembly.
///
/// The compiler will place the branch that is the 
/// most likely just branching instruction that depend
/// on the condition
/// 
/// To look at the assembly:
///
/// ```sh
/// cargo rustc --release --lib --features check_assembly -- --emit=asm
/// cat target/release/deps/*.s | c++filt > readme.S
/// vim readme.S
///  ```
///
/// # Conclusion
///
/// - Macros and likely unlikely functions provided by this crate works
/// - `#[cold]` attribute may inform the compiler about the likelihood of
/// a branch if and only if the function is not inlined by the compiler. To
/// ensure that the cold attribute actualy works, the function must be attributed with
/// `#[inline(never)]`.
pub mod check_assembly {
    extern "C" {
        fn func1(_:u32,_:u32) -> u32;
        fn func2(_:u32,_:u32) -> u32;

        #[cold]
        fn ext_func_cold(_:u32,_:u32) -> u32;
        fn ext_func_hot(_:u32,_:u32) -> u32;
    }

    #[cold]
    fn inlinable_func_cold(a:u32,b:u32) -> u32 {
       unsafe{func1(a,b) }
    }
    fn inlinable_func_hot_a(a:u32,b:u32) -> u32 {
       unsafe{func2(a,b) }
    }

    #[cold]
    #[inline(never)]
    fn inline_never_func_cold(a:u32,b:u32) -> u32 {
       unsafe{func1(a,b) }
    }
    #[inline(never)]
    fn inline_never_func_hot(a:u32,b:u32) -> u32 {
       unsafe{func2(a,b) }
    }
   
    #[inline(never)]
    pub fn check_if_likely_1(a:u32,b:u32,c:u32) -> u32 {
        if_likely!( a<b =>  {
            unsafe{func1(a,b)}//non jumping branch => Ok
        } else {
            unsafe{func2(b,c)}
        })
    }
    #[inline(never)]
    pub fn check_if_likely_2(a:u32,b:u32,c:u32) -> u32 {
        if_likely!( a<b =>  {
            unsafe{func2(b,c)}//non jumping branch => Ok
        } else {
            unsafe{func1(a,b)}
        })
    }
    #[inline(never)]
    pub fn check_if_unlikely_2(a:u32, b:u32,c:u32) -> u32 {
        if_unlikely!( a<b =>  {
            unsafe{func1(a,b)}
        } else {
            unsafe{func2(b,c)}//non jumping branch => Ok
        })
    }
    #[inline(never)]
    pub fn check_if_unlikely_1(a:u32, b:u32,c:u32) -> u32 {
        if_unlikely!( a<b =>  {
            unsafe{func2(b,c)}
        } else {
            unsafe{func1(a,b)}//non jumping branch => Ok
        })
    }
    #[inline(never)]
    pub fn check_if_likely_let_1(a:Option<u32>,b:u32,c:u32) -> u32 {
        if_likely!( let Some(_a) = a =>  {
            unsafe{func1(_a,b)}//non jumping branch => Ok
        } else {
            unsafe{func2(b,c)}
        })
    }
    #[inline(never)]
    pub fn check_if_likely_let_2(a:Option<u32>,b:u32,c:u32) -> u32 {
        if_likely!( let Some(_a) = a =>  {
            unsafe{func2(_a,b)}//non jumping branch => Ok
        } else {
            unsafe{func1(b,c)}
        })
    }
    #[inline(never)]
    pub fn check_if_unlikely_let_2(a:Option<u32>, b:u32,c:u32) -> u32 {
        if_unlikely!( let Some(_a) = a =>  {
            unsafe{func1(_a,b)}
        } else {
            unsafe{func2(b,c)}//non jumping branch => Ok
        })
    }
    #[inline(never)]
    pub fn check_if_unlikely_let_1(a:Option<u32>, b:u32,c:u32) -> u32 {
        if_unlikely!( let Some(_a) = a =>  {
            unsafe{func2(_a,b)}
        } else {
            unsafe{func1(b,c)}//non jumping branch => Ok
        })
    }

    //Here it looks like cold may work ...
    #[inline(never)]
    pub fn check_ext_func_cold_1(a:u32,b:u32,c:u32) ->u32 {
        if a<b  {
            unsafe{ext_func_cold(a,b)}
        } else {
            unsafe{ext_func_hot(b,c)}//non jumping branch => Ok
        }
    }
    //Here it looks like cold may work ...
    #[inline(never)]
    pub fn check_ext_func_cold_2(a:u32,b:u32,c:u32) ->u32 {
        if a<b  {
            unsafe{ext_func_hot(b,c)}//non jumping branch => Ok
        } else {
            unsafe{ext_func_cold(a,b)}
        }
    }

    // but coldness is not preserved if the fonctions is inlined
    #[inline(never)]
    pub fn check_inlinable_func_cold_1(a:u32,b:u32,c:u32) ->u32 {
        if a<b  {
            inlinable_func_cold(a,b)//non jumping branch => NOT OK
        } else {
            inlinable_func_hot_a(b,c)
        }
    }
    #[inline(never)]
    pub fn check_inlinable_func_cold_2(a:u32,b:u32,c:u32) -> u32{
        if a<b  {
            inlinable_func_hot_a(b,c)//non jumping branch => Ok
        } else {
            inlinable_func_cold(a,b)
        }
    }

    // coldness only work for non inlined functions
    #[inline(never)]
    pub fn check_non_inlinable_func_cold_1(a:u32,b:u32,c:u32) ->u32 {
        if a<b  {
            inline_never_func_cold(a,b)
        } else {
            inline_never_func_hot(b,c)//non jumping branch => Ok
        }
    }
    #[inline(never)]
    pub fn check_non_inlinable_func_cold_2(a:u32,b:u32,c:u32) -> u32{
        if a<b  {
            inline_never_func_hot(b,c)//non jumping branch => Ok
        } else {
            inline_never_func_cold(a,b)
        }
    }
}

use std::iter::FromIterator;
use std::mem::MaybeUninit;

pub enum CollectArrayOutput<T, const N: usize> {
    Ok([T; N]),
    Err,
}

impl<T, const N: usize> FromIterator<T> for CollectArrayOutput<T, N> {
    fn from_iter<U: IntoIterator<Item = T>>(iter: U) -> Self {
        let mut iter = iter.into_iter();

        struct Guard<T, const N: usize> {
            ptr: *mut T,
            cnt: usize,
        }

        impl<T, const N: usize> Drop for Guard<T, N> {
            fn drop(&mut self) {
                let initialized_part = std::ptr::slice_from_raw_parts_mut(self.ptr, self.cnt);
                unsafe {
                    std::ptr::drop_in_place(initialized_part);
                }
            }
        }

        let mut arr: [MaybeUninit<T>; N] = unsafe { MaybeUninit::uninit().assume_init() };
        let mut guard = Guard::<T, N> {
            ptr: arr.as_mut_ptr() as *mut T,
            cnt: 0,
        };

        for el in &mut arr {
            match iter.next() {
                Some(v) => {
                    el.write(v);
                    guard.cnt += 1;
                }
                None => {
                    return Self::Err;
                }
            }
        }

        if iter.next().is_some() {
            return Self::Err;
        }

        std::mem::forget(guard);
        Self::Ok(arr.map(|x| unsafe { x.assume_init() }))
    }
}

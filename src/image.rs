use std::mem;

pub struct ImageF {
    pub xsize_: usize,
    pub ysize_: usize,
    pub bytes_per_row: usize,
}

impl ImageF {
    fn bytes_per_row(xsize: usize) -> usize {
        let row_size = xsize + mem::size_of::<f64>() + 32;

        let align = 64;

        let mut bytes_per_row = (row_size + align - 1) & !(align - 1);

        if bytes_per_row % 2048 == 0 {
            bytes_per_row += align;
        }

        bytes_per_row
    }
}

// add others from butteraugli.rs
